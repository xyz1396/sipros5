#include "cfgParser.h"
#include <cmath>
#include <iomanip>
#include <set>

string configPath = "";
string outPath = "";
string element = "";
set<string> supportedElements = {"C13", "H2", "N15", "O18", "S33", "S34"};
double sipLowerBound = 0, sipUpperBound = 100, sipStep = 1;
string help = "Usage: \n \
    -i config file input path\n \
    -o config files output path\n \
    -e SIP label element, e.g., C13, H2, N15, O18, S33, S34\n \
    -l lower bound of SIP abundance\n \
    -u upper bound of SIP abundance\n \
    -s step of SIP abundance \n \
    -h or --help for help";

bool parseArgs(int argc, char const *argv[])
{
    vector<string> vsArguments;
    while (argc--)
        vsArguments.push_back(*argv++);
    for (int i = 1; i < (int)vsArguments.size(); i++)
    {
        if (vsArguments[i] == "-i")
        {
            i = i + 1;
            if (i < (int)vsArguments.size())
            {
                configPath = vsArguments[i];
            }
            else
            {
                cout << help << endl;
                return false;
            }
        }
        else if (vsArguments[i] == "-o")
        {
            i = i + 1;
            if (i < (int)vsArguments.size())
            {
                outPath = vsArguments[i];
            }
            else
            {
                cout << help << endl;
                return false;
            }
        }
        else if (vsArguments[i] == "-e")
        {
            i = i + 1;
            if (i < (int)vsArguments.size())
            {
                element = vsArguments[i];
            }
            else
            {
                cout << help << endl;
                return false;
            }
            if (supportedElements.find(element) == supportedElements.end())
            {
                cout << "Elemement: " << element << " not support" << endl;
                cout << help << endl;
                return false;
            }
        }
        else if (vsArguments[i] == "-l")
        {
            i = i + 1;
            if (i < (int)vsArguments.size())
            {
                sipLowerBound = stof(vsArguments[i]);
            }
            else
            {
                cout << help << endl;
                return false;
            }
        }
        else if (vsArguments[i] == "-u")
        {
            i = i + 1;
            if (i < (int)vsArguments.size())
            {
                sipUpperBound = stof(vsArguments[i]);
            }
            else
            {
                cout << help << endl;
                return false;
            }
        }
        else if (vsArguments[i] == "-s")
        {
            i = i + 1;
            if (i < (int)vsArguments.size())
            {
                sipStep = stof(vsArguments[i]);
            }
            else
            {
                cout << help << endl;
                return false;
            }
        }
        else if (vsArguments[i] == "-h" || vsArguments[i] == "--help")
        {
            cout << help << endl;
            return false;
        }
    }
    if (configPath.length() == 0 || outPath.length() == 0 || element.length() == 0)
    {
        cout << help << endl;
        return false;
    }
    return true;
}

std::string doubleToStringWithWidth(double value, int width, int precision)
{
    std::ostringstream out;
    out << std::fixed << std::setw(width) << std::setfill('0') << std::setprecision(precision) << value;
    return out.str();
}

bool generateCFGs(string cfgPath, string outPath, string element)
{
    cfgParser parser(cfgPath);
    parser.setSearch_NameIX();
    parser.setParent_Mass_WindowsIX();
    vector<int> centers, widths;
    vector<double> pcts;
    for (double pct = sipLowerBound; pct <= sipUpperBound; pct += sipStep)
    {
        pcts.push_back(pct);
        centers.push_back(0);
        // function for width of mass window shift
        int width = round((1.0 - 2.0 * (pct / 100.0 - 0.5) * (pct / 100.0 - 0.5)) * 5.0);
        if (width < 2)
            width = 2;
        widths.push_back(width);
    }
    // in case of upper bound is not included
    pcts.back() = sipUpperBound;
    parser.setElement_PercentIX(element);
    for (size_t i = 0; i < pcts.size(); i++)
    {
        parser.changeSearchName(doubleToStringWithWidth(pcts[i], 6, 3) + "Pct", element);
        parser.changeMassWindowsCenter(centers[i], widths[i]);
        parser.changeSIPabundance(pcts[i], element);
        parser.writeFile(outPath);
    }
    return true;
}

int main(int argc, char const *argv[])
{
    if (parseArgs(argc, argv))
    {
        generateCFGs(configPath, outPath, element);
    }
    return 0;
}
