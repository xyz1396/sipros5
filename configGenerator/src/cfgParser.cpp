#include "cfgParser.h"

cfgParser::cfgParser(const string &cfgFileName)
{
    if (fs::exists(cfgFileName))
    {
        cfgFileStream.open(cfgFileName.c_str(), ios::in);
        if (!cfgFileStream.is_open())
        {
            cout << "Cannot open " << cfgFileName << endl;
        }
        string line;
        while (getline(cfgFileStream, line))
        {
            lines.push_back(line);
        }
        parseParameters();
    }
    else
        cout << cfgFileName << " does not exists" << endl;
}

cfgParser::~cfgParser()
{
    if (cfgFileStream.is_open())
        cfgFileStream.close();
}

void cfgParser::splitString(const string &mString)
{
    tokens.clear();
    size_t start = 0;
    size_t end = 0;
    while (end < mString.size())
    {
        if (mString[end] == '\t' || mString[end] == ' ')
        {
            // ignore continuous sep
            if (start < end)
                tokens.push_back(mString.substr(start, end - start));
            start = end + 1;
        }
        end++;
    }
    // save the last token
    tokens.push_back(mString.substr(start));
}

void cfgParser::parseParameters()
{
    for (size_t i = 0; i < lines.size(); i++)
    {
        if (lines[i][0] != '#')
        {
            splitString(lines[i]);
            parametersIXMap.insert({tokens[0], i});
            parametersMap.insert({tokens[0], tokens});
        }
    }
}

size_t cfgParser::findParameter(const string &parameter)
{
    auto iter = parametersIXMap.find(parameter);
    if (iter != parametersIXMap.end())
        return iter->second;
    else
    {
        return 0;
        cout << parameter << "Not Found!" << endl;
    }
}

void cfgParser::setSearch_NameIX()
{
    Search_NameIX = findParameter("Search_Name");
}

void cfgParser::setParent_Mass_WindowsIX()
{
    Parent_Mass_WindowsIX = findParameter("Parent_Mass_Windows");
}

void cfgParser::setElement_PercentIX(const string &element)
{
    // in case of element is "S33", "S34"
    string token = "Element_Percent{" + element.substr(0, 1) + "}";
    Element_PercentName = token;
    Element_PercentIX = findParameter(token);
}

void cfgParser::changeSearchName(const string &nameSuffix, const string &element)
{
    newFileName = parametersMap["Search_Name"][2] + "_" + element + "_" + nameSuffix;
    lines[Search_NameIX] = "Search_Name = " + newFileName;
    size_t SIP_Element_IsotopeIX = findParameter("SIP_Element_Isotope");
    lines[SIP_Element_IsotopeIX] = "SIP_Element_Isotope = " + element.substr(1);
    size_t SIP_ElementIX = findParameter("SIP_Element");
    lines[SIP_ElementIX] = "SIP_Element = " + element.substr(0, 1);
}

void cfgParser::changeMassWindowsCenter(const int center, const int windowsSize)
{
    vector<int> leftWindows, rightWindows;
    for (int i = 0; i <= windowsSize; i++)
    {
        leftWindows.push_back(center - i);
        rightWindows.push_back(center + i);
    }
    reverse(leftWindows.begin(), leftWindows.end());
    leftWindows.insert(leftWindows.end(), rightWindows.begin() + 1, rightWindows.end());
    string windows;
    for (size_t i = 0; i < leftWindows.size(); i++)
    {
        windows += to_string(leftWindows[i]) + ",";
    }
    windows.pop_back();
    lines[Parent_Mass_WindowsIX] = "Parent_Mass_Windows = " + windows;
}

void cfgParser::changeSIPabundance(const double sipAbundance, const string &element)
{
    string pct;
    double pctD = 0;
    string labelPCT = to_string_with_precision(sipAbundance / 100.0, 5);
    if (element == "O18")
    {
        pctD = (100.0 - sipAbundance) / 100.0 - 0.00038;
        if (pctD < 0)
        {
            cout << "SIP abundance is maxium 99.9 for O18" << endl;
            exit(0);
        }
        pct = to_string_with_precision(pctD, 5);
        pct += ",\t0.00038,\t" + labelPCT;
    }
    else if (element == "S33")
    {
        pctD = (100.0 - sipAbundance) / 100.0 - 0.0429 - 0.0000 - 0.0002;
        if (pctD < 0)
        {
            cout << "SIP abundance is maxium 95.5 for S33" << endl;
            exit(0);
        }
        pct = to_string_with_precision(pctD, 5);
        pct += ",\t" + labelPCT + ",\t0.0429,\t0.0000,\t0.0002";
    }
    else if (element == "S34")
    {
        pctD = (100.0 - sipAbundance) / 100.0 - 0.0076 - 0.0000 - 0.0002;
        if (pctD < 0)
        {
            cout << "SIP abundance is maxium 99.2 for S34" << endl;
            exit(0);
        }
        pct = to_string_with_precision(pctD, 5);
        pct += ",\t0.0076,\t" + labelPCT + ",\t0.0000,\t0.0002";
    }
    else
    {
        pct = to_string_with_precision((100.0F - sipAbundance) / 100.0, 5);
        pct += ",\t" + labelPCT;
    }
    lines[Element_PercentIX] = Element_PercentName + " \t=\t" + pct;
}

void cfgParser::writeFile(const string &folderPath)
{
    fs::path path{folderPath};
    if (!fs::exists(path))
    {
        cout << path.string() << "Not exists" << endl;
        cout << "Creat" << path.string() << endl;
        fs::create_directories(path);
    }
    path /= newFileName + ".cfg";
    std::ofstream out(path);
    for (string &line : lines)
    {
        out << line << "\n";
    }
    out.close();
}