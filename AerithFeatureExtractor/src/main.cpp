#include "PSMfeatureExtractor.h"
#include <omp.h>
#include <getopt.h>

void extractPSMfeaturesTargetAndDecoytoPercolatorPin(const std::string &targetPath,
                                                     const std::string &decoyPath,
                                                     const int topN,
                                                     const std::string &ftFilepath,
                                                     const std::string &configPath,
                                                     const int ThreadNumber,
                                                     const bool doProteinInference,
                                                     const std::string &fileName)
{
    ProNovoConfig::setFilename(configPath);
    PSMfeatureExtractor extractor;
    extractor.extractPSMfeatureParallel(targetPath, decoyPath, topN, ftFilepath, ThreadNumber);
    extractor.writePecorlatorPin(fileName, doProteinInference);
}

std::string helpInfo = R"(
Convert target and decoy PSMs to Percolator pin file 
Usage:
-t <targetPath> -d <decoyPath> -n <topN> 
-f <ftFilepath> -j <ThreadNumber> -c <configPath>
-p <doProteinInference> -o <outputFilePrefix>
)";

int main(int argc, char *argv[])
{
    std::string targetPath, decoyPath, ftFilepath, fileName = "a.pin";
    int topN = 5, ThreadNumber = 3;
    bool doProteinInference = false;
    std::string configPath = "SIP.cfg";
    int opt;
    while ((opt = getopt(argc, argv, "t:d:n:f:c:j:p:o:")) != -1)
    {
        switch (opt)
        {
        case 't':
            targetPath = optarg;
            break;
        case 'd':
            decoyPath = optarg;
            break;
        case 'n':
            topN = std::stoi(optarg);
            break;
        case 'f':
            ftFilepath = optarg;
            break;
        case 'c':
            configPath = optarg;
            break;
        case 'j':
            ThreadNumber = std::stoi(optarg);
            break;
        case 'p':
            doProteinInference = std::stoi(optarg);
            break;
        case 'o':
            fileName = optarg;
            break;
        default:
            std::cerr << helpInfo << std::endl;
            return 1;
        }
    }
    if (targetPath.empty() || decoyPath.empty() || topN == 0 || ftFilepath.empty())
    {
        std::cerr << helpInfo << std::endl;
        return 1;
    }
    extractPSMfeaturesTargetAndDecoytoPercolatorPin(targetPath, decoyPath, topN,
                                                    ftFilepath, configPath, ThreadNumber,
                                                    doProteinInference, fileName);
    return 0;
}