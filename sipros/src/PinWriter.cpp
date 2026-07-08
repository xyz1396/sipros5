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

void PinWriter::writePecorlatorPin(const std::string &fileName, const std::vector<sipPSM> &psms, bool doProteinInference)
{
    setlocale(LC_ALL, "C");
    std::ios_base::sync_with_stdio(false);
    std::ofstream file(fileName);
    if (!file)
    {
        std::cerr << "Unable to open file for writing.\n";
    }
    file << std::fixed << std::setprecision(6);

    const size_t chunkSize = 10000;
    std::stringstream ss;
    ss << "SpecId"
       << "\t"
       << "Label"
       << "\t"
       << "ScanNr"
       << "\t"
       << "ExpMass"
       << "\t"
       << "retentiontime"
       << "\t"
       << "ranks"
       << "\t"
       << "parentCharges"
       << "\t"
       << "massErrors"
       << "\t"
       << "isotopicMassWindowShifts"
       << "\t"
       << "mzShiftFromisolationWindowCenters"
       << "\t"
       << "peptideLengths"
       << "\t"
       << "missCleavageSiteNumbers"
       << "\t"
       << "PTMnumbers"
       << "\t"
       << "isotopicPeakNumbers"
       << "\t"
       << "MS1IsotopicAbundances"
       << "\t"
       << "MS2IsotopicAbundances"
       << "\t"
       << "isotopicAbundanceDiffs"
       << "\t"
       << "WDPscores"
       << "\t"
       << "XcorrScores"
       << "\t"
       << "MVHscores"
       << "\t"
       << "diffScores"
       << "\t"
       << "log10_precursorIntensities"
       << "\t"
       << "Peptide"
       << "\t"
       << "Proteins"
       << "\n";
    file << ss.str();
    ss.str(std::string());
    ss.clear();
    std::string proteinName;
    std::string peptideSeq;
    for (size_t i = 0; i < psms.size(); i++)
    {
        for (size_t j = 0; j < psms[i].fileNames.size(); j += chunkSize)
        {
            for (size_t k = j; k < std::min(j + chunkSize, psms[i].fileNames.size()); ++k)
            {
                ss << psms[i].fileNames[k] << "." << psms[i].scanNumbers[k] << "." << psms[i].ranks[k] << "\t";
                ss << (psms[i].isDecoys[k] ? -1 : 1) << "\t";
                ss << psms[i].scanNumbers[k] << "\t"
                   << psms[i].calculatedParentMasses[k] << "\t"
                   << psms[i].retentionTimes[k] << "\t"
                   << psms[i].ranks[k] << "\t"
                   << psms[i].parentCharges[k] << "\t"
                   << psms[i].massErrors[k] << "\t"
                   << psms[i].isotopicMassWindowShifts[k] << "\t"
                   << psms[i].mzShiftFromisolationWindowCenters[k] << "\t"
                   << psms[i].peptideLengths[k] << "\t"
                   << psms[i].missCleavageSiteNumbers[k] << "\t"
                   << psms[i].PTMnumbers[k] << "\t"
                   << psms[i].isotopicPeakNumbers[k] << "\t"
                   << psms[i].MS1IsotopicAbundances[k] << "\t"
                   << psms[i].MS2IsotopicAbundances[k] << "\t"
                   << psms[i].isotopicAbundanceDiffs[k] << "\t"
                   << psms[i].WDPscores[k] << "\t"
                   << psms[i].XcorrScores[k] << "\t"
                   << psms[i].MVHscores[k] << "\t"
                   << psms[i].diffScores[k] << "\t"
                   << (psms[i].precursorIntensities[k] > 0
                           ? std::log10(psms[i].precursorIntensities[k])
                           : 0)
                   << "\t";

                if (doProteinInference)
                {
                    // for percolator protein inference format
                    peptideSeq = psms[i].originalPeptides[k];
                    if (peptideSeq.front() == '[')
                        peptideSeq.insert(peptideSeq.begin(), 'n');
                    if (peptideSeq.back() == ']')
                        peptideSeq.push_back('n');
                    std::replace(peptideSeq.begin(), peptideSeq.end(), '[', '.');
                    std::replace(peptideSeq.begin(), peptideSeq.end(), ']', '.');
                }
                else
                    peptideSeq = psms[i].identifiedPeptides[k];
                ss << peptideSeq << "\t";

                proteinName = psms[i].proteinNames[k];
                if (doProteinInference)
                {
                    // for percolator protein inference format
                    proteinName = proteinName.substr(1, psms[i].proteinNames[k].length() - 2);
                    std::replace(proteinName.begin(), proteinName.end(), ',', '\t');
                }
                ss << proteinName;
                ss << "\n";
            }
            file << ss.str();
            ss.str(std::string()); // Clear the stringstream
        }
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
        "MS1IsotopicAbundances\tMS2IsotopicAbundances\tisotopicAbundanceDiffs\t"
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
