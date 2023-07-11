/*
 * SiprosReader.cpp
 *
 *  Created on: May 16, 2017
 *      Author: xgo
 */
// fix mzml reading bug, Yi Xiong, 07/10/2023

#include <string>
#include <iostream>
#include "SiprosReader.h"

bool SiprosReader::MzmlReader(string &filename_str, vector<Spectrum> *_vSpectra)
{
	bool bSearchSucceeded = ReadSpectrumData(filename_str.c_str(), _vSpectra);
	return bSearchSucceeded;
}

static void SetMSLevelFilter(MSReader &mstReader)
{
	vector<MSSpectrumType> msLevel;
	msLevel.push_back(MS2);
	mstReader.setFilter(msLevel);
}

bool SiprosReader::ReadSpectrumData(const char *cFile, vector<Spectrum> *_vSpectra)
{
	bool bSucceeded = true;

	// For file access using MSToolkit.
	MSReader mstReader;

	// We want to read only MS2 scans.
	SetMSLevelFilter(mstReader);

	bSucceeded = LoadSpectra(mstReader, cFile, _vSpectra);

	return bSucceeded;
}

bool SiprosReader::LoadSpectra(MSReader &mstReader, const char *cFile, vector<Spectrum> *_vSpectra)
{
	// For holding spectrum.
	Spectrum mstSpectrum;
	int lastScanNumber, currentScanNumber;
	// read first MS2 scan
	PreloadIonsFile(mstReader, mstSpectrum, cFile, false);
	lastScanNumber = mstReader.getLastScan();
	_vSpectra->reserve(lastScanNumber);
	currentScanNumber = mstSpectrum.getScanNumber();
	// Load all MS2 scans. If reach over last scan, scan number will be 0
	while (currentScanNumber != 0 && currentScanNumber <= lastScanNumber)
	{
		_vSpectra->push_back(mstSpectrum);
		// read next scan
		PreloadIonsFile(mstReader, mstSpectrum, NULL, true);
		currentScanNumber = mstSpectrum.getScanNumber();
	}
	return true;
}

void SiprosReader::PreloadIonsFile(MSReader &mstReader, Spectrum &spec, const char *cFile, bool bNext, int scNum)
{
	// open file for first scan
	if (!bNext)
		mstReader.readFile(cFile, spec, scNum);
	// for next scans
	else
		mstReader.readFile(NULL, spec);
}
