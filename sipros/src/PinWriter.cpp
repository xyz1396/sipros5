#include "PinWriter.h"

#include <algorithm>
#include <clocale>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <locale>
#include <map>
#include <sstream>
#include <utility>

#include <omp.h>

namespace
{

constexpr size_t PinChunkSize = 10000;

struct PinFormatChunk
{
    size_t psmIndex = 0;
    size_t begin = 0;
    size_t end = 0;
    std::string text;
};

std::string formatPsmChunk(const sipPSM &psm,
                           size_t begin,
                           size_t end,
                           bool doProteinInference)
{
    std::ostringstream ss;
    std::string proteinName;
    std::string peptideSeq;
    for (size_t k = begin; k < end; ++k)
    {
        ss << psm.fileNames[k] << "." << psm.scanNumbers[k] << "." << psm.ranks[k] << "\t";
        ss << (psm.isDecoys[k] ? -1 : 1) << "\t";
        ss << psm.scanNumbers[k] << "\t"
           << psm.calculatedParentMasses[k] << "\t"
           << psm.retentionTimes[k] << "\t"
           << psm.ranks[k] << "\t"
           << psm.parentCharges[k] << "\t"
           << psm.massErrors[k] << "\t"
           << psm.isotopicMassWindowShifts[k] << "\t"
           << psm.mzShiftFromisolationWindowCenters[k] << "\t"
           << psm.peptideLengths[k] << "\t"
           << psm.missCleavageSiteNumbers[k] << "\t"
           << psm.PTMnumbers[k] << "\t"
           << psm.isotopicPeakNumbers[k] << "\t"
           << psm.MS1IsotopeFitScores[k] << "\t"
           << psm.MS1IsotopicAbundances[k] << "\t"
           << psm.MS2IsotopicAbundances[k] << "\t"
           << psm.isotopicAbundanceDiffs[k] << "\t"
           << psm.WDPscores[k] << "\t"
           << psm.XcorrScores[k] << "\t"
           << psm.MVHscores[k] << "\t"
           << psm.diffScores[k] << "\t"
           << psm.matchedBIons[k] << "\t"
           << psm.matchedYIons[k] << "\t"
           << psm.maxConsecutiveBIons[k] << "\t"
           << psm.maxConsecutiveYIons[k] << "\t"
           << (psm.precursorIntensities[k] > 0
                   ? std::log10(psm.precursorIntensities[k])
                   : 0)
           << "\t";

        if (doProteinInference)
        {
            peptideSeq = psm.originalPeptides[k];
            if (peptideSeq.front() == '[')
                peptideSeq.insert(peptideSeq.begin(), 'n');
            if (peptideSeq.back() == ']')
                peptideSeq.push_back('n');
            std::replace(peptideSeq.begin(), peptideSeq.end(), '[', '.');
            std::replace(peptideSeq.begin(), peptideSeq.end(), ']', '.');
        }
        else
        {
            peptideSeq = psm.identifiedPeptides[k];
        }
        ss << peptideSeq << "\t";

        proteinName = psm.proteinNames[k];
        if (doProteinInference)
        {
            proteinName = proteinName.substr(1, psm.proteinNames[k].length() - 2);
            std::replace(proteinName.begin(), proteinName.end(), ',', '\t');
        }
        ss << proteinName << "\n";
    }
    return ss.str();
}

} // namespace

void PinWriter::writePecorlatorPin(const std::string &fileName, const std::vector<sipPSM> &psms, bool doProteinInference)
{
    setlocale(LC_ALL, "C");
    std::ios_base::sync_with_stdio(false);
    static const char *header =
        "SpecId\tLabel\tScanNr\tExpMass\tretentiontime\tranks\tparentCharges\t"
        "massErrors\tisotopicMassWindowShifts\tmzShiftFromisolationWindowCenters\t"
        "peptideLengths\tmissCleavageSiteNumbers\tPTMnumbers\tisotopicPeakNumbers\t"
        "MS1IsotopeFitScore\tMS1IsotopicAbundances\tMS2IsotopicAbundances\t"
        "isotopicAbundanceDiffs\tWDPscores\tXcorrScores\tMVHscores\tdiffScores\t"
        "matchedBIons\tmatchedYIons\tmaxConsecutiveBIons\tmaxConsecutiveYIons\t"
        "log10_precursorIntensities\tPeptide\tProteins\n";

    std::vector<PinFormatChunk> chunks;
    for (size_t i = 0; i < psms.size(); ++i)
    {
        for (size_t begin = 0; begin < psms[i].fileNames.size(); begin += PinChunkSize)
        {
            chunks.push_back(
                {i, begin, std::min(begin + PinChunkSize, psms[i].fileNames.size()), {}});
        }
    }

#pragma omp parallel for schedule(static)
    for (int64_t i = 0; i < static_cast<int64_t>(chunks.size()); ++i)
    {
        PinFormatChunk &chunk = chunks[static_cast<size_t>(i)];
        chunk.text = formatPsmChunk(
            psms[chunk.psmIndex], chunk.begin, chunk.end, doProteinInference);
    }

    std::ofstream file(fileName, std::ios::binary);
    if (!file)
    {
        std::cerr << "Unable to open file for writing.\n";
        return;
    }
    file.write(header, static_cast<std::streamsize>(std::char_traits<char>::length(header)));
    for (const PinFormatChunk &chunk : chunks)
    {
        file.write(chunk.text.data(), static_cast<std::streamsize>(chunk.text.size()));
    }
    file.close();
}


size_t PinWriter::writeSearchSpectraPin(const std::string &fileName,
                                        const std::string &sampleBasename,
                                        const std::vector<SearchSpectraPinRow> &rows)
{
    static const char *header =
        "SpecId\tLabel\tScanNr\tExpMass\tretentiontime\tranks\tparentCharges\t"
        "massErrors\tisotopicMassWindowShifts\tmzShiftFromisolationWindowCenters\t"
        "peptideLengths\tmissCleavageSiteNumbers\tPTMnumbers\tisotopicPeakNumbers\t"
        "MS1IsotopeFitScore\tMS1IsotopicAbundances\tMS2IsotopicAbundances\t"
        "isotopicAbundanceDiffs\t"
        "WDPscores\tXcorrScores\tMVHscores\tentropyScores\tcosineScores\t"
        "matchedYenvelopes\tmatchedBenvelopes\tdeltaRT\t"
        "diffScores\tlog10_precursorIntensities\tPeptide\tProteins\n";

    std::map<int, double> topByScan;
    for (const SearchSpectraPinRow &row : rows)
    {
        auto it = topByScan.find(row.scanNumber);
        if (it == topByScan.end() || row.wdpScore > it->second)
        {
            topByScan[row.scanNumber] = row.wdpScore;
        }
    }

    std::ostringstream pin;
    pin.imbue(std::locale::classic());
    pin << header;
    pin << std::fixed << std::setprecision(6);
    for (const SearchSpectraPinRow &row : rows)
    {
        const double top = topByScan[row.scanNumber];
        const double diff = top - row.wdpScore;
        pin << sampleBasename << "." << row.scanNumber << "." << row.rank << '\t';
        pin << row.label << '\t';
        pin << row.scanNumber << '\t';
        pin << row.expMass << '\t';
        pin << row.retentionTime << '\t';
        pin << row.rank << '\t';
        pin << row.parentCharge << '\t';
        pin << row.massError << '\t';
        pin << row.isotopicMassWindowShift << '\t';
        pin << row.mzShiftFromIsolationWindowCenter << '\t';
        pin << row.peptideLength << '\t';
        pin << row.missCleavageSiteNumber << '\t';
        pin << row.ptmNumber << '\t';
        pin << row.isotopicPeakNumber << '\t';
        pin << row.ms1IsotopeFitScore << '\t';
        pin << row.ms1IsotopicAbundance << '\t';
        pin << row.ms2IsotopicAbundance << '\t';
        pin << (row.ms1IsotopicAbundance - row.ms2IsotopicAbundance) << '\t';
        pin << row.wdpScore << '\t';
        pin << row.xcorrScore << '\t';
        pin << row.mvhScore << '\t';
        pin << row.entropyScore << '\t';
        pin << row.cosineScore << '\t';
        pin << row.matchedYEnvelope << '\t';
        pin << row.matchedBEnvelope << '\t';
        pin << row.deltaRT << '\t';
        pin << diff << '\t';
        pin << row.log10PrecursorIntensity << '\t';
        pin << row.peptide << '\t';
        pin << row.proteins << '\n';
    }

    const std::string pinText = pin.str();
    std::ofstream os(fileName);
    if (!os)
    {
        std::cerr << "Cannot write " << fileName << "\n";
        return 0;
    }
    os.write(pinText.data(), static_cast<std::streamsize>(pinText.size()));
    if (!os)
    {
        std::cerr << "Failed while writing " << fileName << "\n";
        return 0;
    }
    return rows.size();
}
