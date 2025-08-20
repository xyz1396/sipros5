#include "cfgParser.h"
#include <cmath>
#include <iomanip>
#include <set>

string configPath = "";
string outPath = "";
string element = "";
float toleranceMS1 = 0.01;
float toleranceMS2 = 0.01;
set<string> supportedElements = {"C13", "H2", "N15", "O18", "S33", "S34", "R"};
double sipLowerBound = 0, sipUpperBound = 100, sipStep = 1;
string help = "Usage:\n"
              "  -i <input_config_path>    Path to input config file (required)\n"
              "  -o <output_path>          Path to output config files (required)\n"
              "  -e <element>              SIP label element: C13, H2, N15, O18, S33, S34, R (R = regular search) (required)\n"
              "  -l <lower_bound>          Lower bound of SIP abundance (default: 0)\n"
              "  -u <upper_bound>          Upper bound of SIP abundance (default: 100)\n"
              "  -s <step>                 Step size for SIP abundance (default: 1)\n"
              "  -t1 <ms1_tolerance>       MS1 tolerance in Da (default: 0.01)\n"
              "  -t2 <ms2_tolerance>       MS2 tolerance in Da (default: 0.01)\n"
              "  -h, --help                Show this help message\n";

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
        else if (vsArguments[i] == "-t1")
        {
            i = i + 1;
            if (i < (int)vsArguments.size())
            {
                toleranceMS1 = stof(vsArguments[i]);
            }
            else
            {
                cout << help << endl;
                return false;
            }
        }
        else if (vsArguments[i] == "-t2")
        {
            i = i + 1;
            if (i < (int)vsArguments.size())
            {
                toleranceMS2 = stof(vsArguments[i]);
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

bool generateRegularCFG(string cfgPath, string outPath)
{
    cfgParser parser(cfgPath);
    parser.newFileName = "Regular";
    parser.changeMStolerance(toleranceMS1, toleranceMS2);
    parser.writeFile(outPath);
    return true;
}

bool generateCFGs(string cfgPath, string outPath, string element)
{
    cfgParser parser(cfgPath);
    parser.setSearch_NameIX();
    parser.changeMStolerance(toleranceMS1, toleranceMS2);
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
        parser.changeSearchName(doubleToStringWithWidth(pcts[i], 7, 3) + "Pct", element);
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
        if (element == "R")
        {
            generateRegularCFG(configPath, outPath);
        }
        else
            generateCFGs(configPath, outPath, element);
    }
    return 0;
}
