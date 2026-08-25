/*
 * CometSearch.cpp
 *
 *  Created on: May 13, 2016
 *      Author: xgo
 */

#include "CometSearchMod.h"

int CometSearchMod::iArraySizePreprocess = 0;
int CometSearchMod::iArraySizeScore = 0;
int CometSearchMod::iMaxPercusorCharge = 0;
int CometSearchMod::iDimesion2 = 0;
int CometSearchMod::iMAX_PEPTIDE_LEN = 0;
double CometSearchMod::ProbabilityCutOff = 0.01;

CometSearchMod::CometSearchMod() {

}

CometSearchMod::~CometSearchMod() {

}

/*
 * variables needed from MS2Scan:
 * isMS2HighRes
 * dParentMass
 * iParentChargeState
 * vdMZ
 */
bool CometSearchMod::Preprocess(struct Query *pScoring, MS2Scan * mstSpectrum, double *pdTmpRawData, double *pdTmpFastXcorrData, double *pdTmpCorrelationData,
		double *pdTmpSmoothedSpectrum, double *pdTmpPeakExtracted,
		double maximumArrayMz, bool buildSparseFastXcorrDirectly) {
	int i;
	int x;
	int y;
	// struct msdata pTmpSpData[NUM_SP_IONS];
	struct PreprocessStruct pPre;
	double dInverseBinWidth = 0; // dFragmentBinSize = 0;
	int iMinus17 = 0, iMinus18 = 0;
	//mstSpectrum->isMS2HighRes = false;
	if (mstSpectrum->isMS2HighRes) {
		dInverseBinWidth = ProNovoConfig::dHighResInverseBinWidth;
		iMinus17 = static_cast<int>(ProNovoConfig::precalcMasses.iMinus17HighRes);
		iMinus18 = static_cast<int>(ProNovoConfig::precalcMasses.iMinus18HighRes);
	} else {
		dInverseBinWidth = ProNovoConfig::dLowResInverseBinWidth;
		iMinus17 = static_cast<int>(ProNovoConfig::precalcMasses.iMinus17LowRes);
		iMinus18 = static_cast<int>(ProNovoConfig::precalcMasses.iMinus18LowRes);
	}
	pPre.iHighestIon = 0;
	pPre.dHighestIntensity = 0;

	//MH: Find appropriately sized array cushion based on user parameters. Fixes error found by Patrick Pedrioli for
	// very wide mass tolerance searches (i.e. 500 Da).
	double dCushion = 3.0;
	double arrayUpperMz = mstSpectrum->dParentMass;
	if (maximumArrayMz > 0.0 && std::isfinite(maximumArrayMz)) {
		// The XCorr correlation spectrum is identically zero above the highest
		// observed fragment plus the background-correlation tail. Window-only
		// precursor hypotheses can otherwise size this array to the z=4 neutral
		// mass even though no MS2 data exist in most of that range. Express the
		// tail in bins so this remains safe for both high- and low-resolution MS2.
		const double correlationTailMz =
			static_cast<double>(ProNovoConfig::iXcorrProcessingOffset + 12) /
			dInverseBinWidth;
		arrayUpperMz = std::min(
			arrayUpperMz, maximumArrayMz + correlationTailMz);
	}
	pScoring->_spectrumInfoInternal.iArraySize = (int) ((arrayUpperMz + dCushion + 2.0) * dInverseBinWidth);
	// seg debug b
	if (pScoring->_spectrumInfoInternal.iArraySize > CometSearchMod::iArraySizePreprocess) {
		cout << "Error 90" << endl;
	}
	// seg debug e
	const int preprocessingCharge = mstSpectrum->bUseReactionChargeForScoring
			? mstSpectrum->iParentChargeState
			: mstSpectrum->iMaxCandidateCharge;
	if (preprocessingCharge <= 1) {
		pScoring->_spectrumInfoInternal.iMaxFragCharge = 1;
	} else {
		pScoring->_spectrumInfoInternal.iMaxFragCharge = preprocessingCharge - 1;
	}
	if (preprocessingCharge > ProNovoConfig::iMaxPercusorCharge) {
		cout << "Error 90" << endl;
		exit(1);
	}

	// Initialize only the portion used by this spectrum.  These buffers are
	// allocated for the largest spectrum in the file, but clearing that global
	// extent for every scan turns XCorr preprocessing into avoidable memory
	// traffic when one high-mass precursor hypothesis inflates the maximum.
	const size_t iTmp = static_cast<size_t>(
		pScoring->_spectrumInfoInternal.iArraySize) * sizeof(double);
	//seg debug b
	if (iTmp > (iArraySizePreprocess * sizeof(double))) {
		cout << "Error 2" << endl;
	}
	//seg debug e
	memset(pdTmpRawData, 0, iTmp);
	memset(pdTmpFastXcorrData, 0, iTmp);
	memset(pdTmpCorrelationData, 0, iTmp);
	memset(pdTmpSmoothedSpectrum, 0, iTmp);
	memset(pdTmpPeakExtracted, 0, iTmp);

	// pdTmpRawData is a binned array holding raw data
	if (!LoadIons(pScoring, pdTmpRawData, mstSpectrum, &pPre)) {
		return false;
	}

	if (!buildSparseFastXcorrDirectly) {
		try {
			pScoring->pfFastXcorrData = new float[pScoring->_spectrumInfoInternal.iArraySize]();
		} catch (std::bad_alloc& ba) {
			const string strErrorMsg = " Error - new(pfFastXcorrData[" +
				to_string(pScoring->_spectrumInfoInternal.iArraySize) +
				"]). bad_alloc: " + ba.what() +
				".\nComet ran out of memory. Adjust spectrum_batch_size to mitigate memory use.\n";
			logerr(strErrorMsg.c_str());
			return false;
		}

		if (ProNovoConfig::ionInformation.bUseNeutralLoss
				&& (ProNovoConfig::ionInformation.iIonVal[ION_SERIES_A] || ProNovoConfig::ionInformation.iIonVal[ION_SERIES_B]
						|| ProNovoConfig::ionInformation.iIonVal[ION_SERIES_Y])) {
			try {
				pScoring->pfFastXcorrDataNL = new float[pScoring->_spectrumInfoInternal.iArraySize]();
			} catch (std::bad_alloc& ba) {
				const string strErrorMsg = " Error - new(pfFastXcorrDataNL[" +
					to_string(pScoring->_spectrumInfoInternal.iArraySize) +
					"]). bad_alloc: " + ba.what() +
					".\nComet ran out of memory. Adjust spectrum_batch_size to mitigate memory use.\n";
				logerr(strErrorMsg.c_str());
				return false;
			}
		}
	}

	// Create data for correlation analysis.
	// pdTmpRawData intensities are normalized to 100; pdTmpCorrelationData is windowed
	MakeCorrData(pdTmpRawData, pdTmpCorrelationData, pScoring, &pPre);

	// Make fast xcorr spectrum.
	double dSum = 0.0;
	int iTmpRange = 2 * ProNovoConfig::iXcorrProcessingOffset + 1;
	double dTmp = 1.0 / (double) (iTmpRange - 1);

	dSum = 0.0;
	for (i = 0; i < ProNovoConfig::iXcorrProcessingOffset; i++) {
		dSum += pdTmpCorrelationData[i];
	}
	for (i = ProNovoConfig::iXcorrProcessingOffset; i < pScoring->_spectrumInfoInternal.iArraySize + ProNovoConfig::iXcorrProcessingOffset; i++) {
		if (i < pScoring->_spectrumInfoInternal.iArraySize) {
			dSum += pdTmpCorrelationData[i];
		}
		if (i >= iTmpRange) {
			dSum -= pdTmpCorrelationData[i - iTmpRange];
		}
		//seg debug b
		if ((i - ProNovoConfig::iXcorrProcessingOffset) >= iArraySizePreprocess) {
			cout << "Error 5" << endl;
		}
		//seg debug e
		pdTmpFastXcorrData[i - ProNovoConfig::iXcorrProcessingOffset] = (dSum - pdTmpCorrelationData[i - ProNovoConfig::iXcorrProcessingOffset]) * dTmp;
	}

	if (buildSparseFastXcorrDirectly) {
		const bool useNeutralLoss =
			ProNovoConfig::ionInformation.bUseNeutralLoss
			&& (ProNovoConfig::ionInformation.iIonVal[ION_SERIES_A]
				|| ProNovoConfig::ionInformation.iIonVal[ION_SERIES_B]
				|| ProNovoConfig::ionInformation.iIonVal[ION_SERIES_Y]);
		pScoring->iFastXcorrData =
			pScoring->_spectrumInfoInternal.iArraySize /
			SPARSE_MATRIX_SIZE + 1;
		try {
			pScoring->ppfSparseFastXcorrData =
				new float*[pScoring->iFastXcorrData]();
		} catch (std::bad_alloc& ba) {
			const string strErrorMsg =
				" Error - new(pScoring->ppfSparseFastXcorrData[" +
				to_string(pScoring->iFastXcorrData) +
				"]). bad_alloc: " + ba.what() +
				".\nComet ran out of memory. Adjust spectrum_batch_size to mitigate memory use.\n";
			logerr(strErrorMsg.c_str());
			return false;
		}
		if (useNeutralLoss) {
			pScoring->iFastXcorrDataNL =
				pScoring->_spectrumInfoInternal.iArraySize /
				SPARSE_MATRIX_SIZE + 1;
			try {
				pScoring->ppfSparseFastXcorrDataNL =
					new float*[pScoring->iFastXcorrDataNL]();
			} catch (std::bad_alloc& ba) {
				const string strErrorMsg =
					" Error - new(pScoring->ppfSparseFastXcorrDataNL[" +
					to_string(pScoring->iFastXcorrDataNL) +
					"]). bad_alloc: " + ba.what() +
					".\nComet ran out of memory. Adjust spectrum_batch_size to mitigate memory use.\n";
				logerr(strErrorMsg.c_str());
				return false;
			}
		}

		for (i = 1; i < pScoring->_spectrumInfoInternal.iArraySize; ++i) {
			float fastXcorr = static_cast<float>(
				pdTmpCorrelationData[i] - pdTmpFastXcorrData[i]);
			if (mstSpectrum->isMS2HighRes) {
				int flank = i - 1;
				fastXcorr += static_cast<float>(
					(pdTmpCorrelationData[flank] -
					 pdTmpFastXcorrData[flank]) * 0.5);
				flank = i + 1;
				if (flank < pScoring->_spectrumInfoInternal.iArraySize) {
					fastXcorr += static_cast<float>(
						(pdTmpCorrelationData[flank] -
						 pdTmpFastXcorrData[flank]) * 0.5);
				}
			}
			if (fastXcorr > FLOAT_ZERO || fastXcorr < -FLOAT_ZERO) {
				x = i / SPARSE_MATRIX_SIZE;
				if (pScoring->ppfSparseFastXcorrData[x] == NULL) {
					try {
						pScoring->ppfSparseFastXcorrData[x] =
							new float[SPARSE_MATRIX_SIZE]();
					} catch (std::bad_alloc& ba) {
						const string strErrorMsg =
							" Error - new sparse fast XCorr block. bad_alloc: " +
							string(ba.what()) +
							".\nComet ran out of memory. Adjust spectrum_batch_size to mitigate memory use.\n";
						logerr(strErrorMsg.c_str());
						return false;
					}
				}
				y = i - (x * SPARSE_MATRIX_SIZE);
				pScoring->ppfSparseFastXcorrData[x][y] = fastXcorr;
			}

			if (useNeutralLoss) {
				float fastXcorrNL = fastXcorr;
				int lossBin = i - iMinus17;
				if (lossBin >= 0) {
					fastXcorrNL += static_cast<float>(
						(pdTmpCorrelationData[lossBin] -
						 pdTmpFastXcorrData[lossBin]) * 0.2);
				}
				lossBin = i - iMinus18;
				if (lossBin >= 0) {
					fastXcorrNL += static_cast<float>(
						(pdTmpCorrelationData[lossBin] -
						 pdTmpFastXcorrData[lossBin]) * 0.2);
				}
				if (fastXcorrNL > FLOAT_ZERO ||
					fastXcorrNL < -FLOAT_ZERO) {
					x = i / SPARSE_MATRIX_SIZE;
					if (pScoring->ppfSparseFastXcorrDataNL[x] == NULL) {
						try {
							pScoring->ppfSparseFastXcorrDataNL[x] =
								new float[SPARSE_MATRIX_SIZE]();
						} catch (std::bad_alloc& ba) {
							const string strErrorMsg =
								" Error - new sparse neutral-loss XCorr block. bad_alloc: " +
								string(ba.what()) +
								".\nComet ran out of memory. Adjust spectrum_batch_size to mitigate memory use.\n";
							logerr(strErrorMsg.c_str());
							return false;
						}
					}
					y = i - (x * SPARSE_MATRIX_SIZE);
					pScoring->ppfSparseFastXcorrDataNL[x][y] =
						fastXcorrNL;
				}
			}
		}
		return true;
	}

	pScoring->pfFastXcorrData[0] = 0.0;
	for (i = 1; i < pScoring->_spectrumInfoInternal.iArraySize; i++) {
		double dTmp = pdTmpCorrelationData[i] - pdTmpFastXcorrData[i];

		pScoring->pfFastXcorrData[i] = (float) dTmp;

		// Add flanking peaks if used if it is high resolution
		if (mstSpectrum->isMS2HighRes) {
			int iTmp;

			iTmp = i - 1;
			pScoring->pfFastXcorrData[i] += (float) ((pdTmpCorrelationData[iTmp] - pdTmpFastXcorrData[iTmp]) * 0.5);

			iTmp = i + 1;
			if (iTmp < pScoring->_spectrumInfoInternal.iArraySize)
				pScoring->pfFastXcorrData[i] += (float) ((pdTmpCorrelationData[iTmp] - pdTmpFastXcorrData[iTmp]) * 0.5);
		}

		// If A, B or Y ions and their neutral loss selected, roll in -17/-18 contributions to pfFastXcorrDataNL
		if (ProNovoConfig::ionInformation.bUseNeutralLoss
				&& (ProNovoConfig::ionInformation.iIonVal[ION_SERIES_A] || ProNovoConfig::ionInformation.iIonVal[ION_SERIES_B]
						|| ProNovoConfig::ionInformation.iIonVal[ION_SERIES_Y])) {
			int iTmp;

			pScoring->pfFastXcorrDataNL[i] = pScoring->pfFastXcorrData[i];

			iTmp = i - iMinus17;
			if (iTmp >= 0) {
				pScoring->pfFastXcorrDataNL[i] += (float) ((pdTmpCorrelationData[iTmp] - pdTmpFastXcorrData[iTmp]) * 0.2);
			}

			iTmp = i - iMinus18;
			if (iTmp >= 0) {
				pScoring->pfFastXcorrDataNL[i] += (float) ((pdTmpCorrelationData[iTmp] - pdTmpFastXcorrData[iTmp]) * 0.2);
			}

		}
	}

	/*ofstream outputFile;
	 string sOutputFile = "pScoring->pfFastXcorrDataNL.txt";
	 outputFile.open(sOutputFile.c_str());
	 for (i = 1; i < pScoring->_spectrumInfoInternal.iArraySize; i++) {
	 if (pScoring->pfFastXcorrDataNL[i] > FLOAT_ZERO || pScoring->pfFastXcorrDataNL[i] < -FLOAT_ZERO) {
	 outputFile << i << "\t" << pScoring->pfFastXcorrDataNL[i] << endl;
	 }
	 }
	 outputFile.close();*/

	// Using sparse matrix which means we free pScoring->pfFastXcorrData, ->pfFastXcorrDataNL here
	// If A, B or Y ions and their neutral loss selected, roll in -17/-18 contributions to pfFastXcorrDataNL.
	if (ProNovoConfig::ionInformation.bUseNeutralLoss
			&& (ProNovoConfig::ionInformation.iIonVal[ION_SERIES_A] || ProNovoConfig::ionInformation.iIonVal[ION_SERIES_B]
					|| ProNovoConfig::ionInformation.iIonVal[ION_SERIES_Y])) {
		pScoring->iFastXcorrDataNL = pScoring->_spectrumInfoInternal.iArraySize / SPARSE_MATRIX_SIZE + 1;

		try {
			pScoring->ppfSparseFastXcorrDataNL = new float*[pScoring->iFastXcorrDataNL]();
		} catch (std::bad_alloc& ba) {
			const string strErrorMsg = " Error - new(pScoring->ppfSparseFastXcorrDataNL[" +
				to_string(pScoring->iFastXcorrDataNL) + "]). bad_alloc: " + ba.what() +
				".\nComet ran out of memory. Adjust spectrum_batch_size to mitigate memory use.\n";
			logerr(strErrorMsg.c_str());
			return false;
		}
		for (i = 0; i < pScoring->iFastXcorrDataNL; i++) {
			pScoring->ppfSparseFastXcorrDataNL[i] = NULL;
		}

		for (i = 1; i < pScoring->_spectrumInfoInternal.iArraySize; i++) {
			if (pScoring->pfFastXcorrDataNL[i] > FLOAT_ZERO || pScoring->pfFastXcorrDataNL[i] < -FLOAT_ZERO) {
				x = i / SPARSE_MATRIX_SIZE;
				// seg debug b
				if (x >= pScoring->iFastXcorrDataNL) {
					cout << "Error 91" << endl;
				}
				// seg debug e
				if (pScoring->ppfSparseFastXcorrDataNL[x] == NULL) {
					try {
						pScoring->ppfSparseFastXcorrDataNL[x] = new float[SPARSE_MATRIX_SIZE]();
					} catch (std::bad_alloc& ba) {
						const string strErrorMsg = " Error - new(pScoring->ppfSparseFastXcorrDataNL[" +
							to_string(x) + "][" + to_string(SPARSE_MATRIX_SIZE) +
							"]). bad_alloc: " + ba.what() +
							".\nComet ran out of memory. Adjust spectrum_batch_size to mitigate memory use.\n";
						logerr(strErrorMsg.c_str());
						return false;
					}
					for (y = 0; y < SPARSE_MATRIX_SIZE; y++) {
						pScoring->ppfSparseFastXcorrDataNL[x][y] = 0;
					}
				}
				y = i - (x * SPARSE_MATRIX_SIZE);
				// seg debug b
				if (y >= SPARSE_MATRIX_SIZE) {
					cout << "Error 92" << endl;
				}
				// seg debug e
				pScoring->ppfSparseFastXcorrDataNL[x][y] = pScoring->pfFastXcorrDataNL[i];
			}
		}

		delete[] pScoring->pfFastXcorrDataNL;
		pScoring->pfFastXcorrDataNL = NULL;

	}

	pScoring->iFastXcorrData = pScoring->_spectrumInfoInternal.iArraySize / SPARSE_MATRIX_SIZE + 1;

	//MH: Fill sparse matrix
	try {
		pScoring->ppfSparseFastXcorrData = new float*[pScoring->iFastXcorrData]();
	} catch (std::bad_alloc& ba) {
		const string strErrorMsg = " Error - new(pScoring->ppfSparseFastXcorrData[" +
			to_string(pScoring->iFastXcorrData) + "]). bad_alloc: " + ba.what() +
			".\nComet ran out of memory. Adjust spectrum_batch_size to mitigate memory use.\n";
		logerr(strErrorMsg.c_str());
		return false;
	}
	for (i = 0; i < pScoring->iFastXcorrData; i++) {
		pScoring->ppfSparseFastXcorrData[i] = NULL;
	}

	for (i = 1; i < pScoring->_spectrumInfoInternal.iArraySize; i++) {
		if (pScoring->pfFastXcorrData[i] > FLOAT_ZERO || pScoring->pfFastXcorrData[i] < -FLOAT_ZERO) {
			x = i / SPARSE_MATRIX_SIZE;
			// seg debug b
			if (x >= pScoring->iFastXcorrData) {
				cout << "Error 94" << endl;
			}
			// seg debug e
			if (pScoring->ppfSparseFastXcorrData[x] == NULL) {
				try {
					pScoring->ppfSparseFastXcorrData[x] = new float[SPARSE_MATRIX_SIZE]();
				} catch (std::bad_alloc& ba) {
					const string strErrorMsg = " Error - new(pScoring->ppfSparseFastXcorrData[" +
						to_string(x) + "][" + to_string(SPARSE_MATRIX_SIZE) +
						"]). bad_alloc: " + ba.what() +
						".\nComet ran out of memory. Adjust spectrum_batch_size to mitigate memory use.\n";
					logerr(strErrorMsg.c_str());
					return false;
				}
				for (y = 0; y < SPARSE_MATRIX_SIZE; y++)
					pScoring->ppfSparseFastXcorrData[x][y] = 0;
			}
			y = i - (x * SPARSE_MATRIX_SIZE);
			// seg debug b
			if (y >= SPARSE_MATRIX_SIZE) {
				cout << "Error 95" << endl;
			}
			// seg debug e
			pScoring->ppfSparseFastXcorrData[x][y] = pScoring->pfFastXcorrData[i];
		}
	}

	delete[] pScoring->pfFastXcorrData;
	pScoring->pfFastXcorrData = NULL;

	return true;
}

//  Reads MSMS data file as ASCII mass/intensity pairs.
bool CometSearchMod::LoadIons(struct Query *pScoring, double *pdTmpRawData, MS2Scan * mstSpectrum, struct PreprocessStruct *pPre) {
	double dIon, dIntensity;
	double dInverseBinWidth = mstSpectrum->isMS2HighRes ? (ProNovoConfig::dHighResInverseBinWidth) : (ProNovoConfig::dLowResInverseBinWidth);
	double dOneMinusBinOffset = mstSpectrum->isMS2HighRes ? (ProNovoConfig::dHighResOneMinusBinOffset) : (ProNovoConfig::dLowResOneMinusBinOffset);
	size_t i = 0;
	pScoring->_spectrumInfoInternal.dTotalIntensity = 0;
	while (true) {
		if (i >= mstSpectrum->vdMZ.size()) {
			break;
		}
		dIon = mstSpectrum->vdMZ.at(i);
		dIntensity = mstSpectrum->vdIntensity.at(i);
		i++;

		pScoring->_spectrumInfoInternal.dTotalIntensity += dIntensity;

		if ((dIntensity >= ProNovoConfig::options.dMinIntensity) && (dIntensity > 0.0)) {
			if (dIon < (mstSpectrum->dParentMass + 50.0)) {
				int iBinIon = BINX(dIon, dInverseBinWidth, dOneMinusBinOffset);
				dIntensity = sqrt(dIntensity);

				if (iBinIon > pPre->iHighestIon) {
					pPre->iHighestIon = iBinIon;
				}

				//seg debug b
				if (iBinIon >= iArraySizePreprocess) {
					cout << "Error 3" << endl;
				}
				//seg debug e

				if ((iBinIon < pScoring->_spectrumInfoInternal.iArraySize) && (dIntensity > pdTmpRawData[iBinIon])) {
					if (ProNovoConfig::options.iRemovePrecursor == 1) {
						double dMZ = (mstSpectrum->dParentMass + (pScoring->_spectrumInfoInternal.iChargeState - 1) * PROTON_MASS)
								/ (double) (pScoring->_spectrumInfoInternal.iChargeState);
						if (fabs(dIon - dMZ) > ProNovoConfig::options.dRemovePrecursorTol) {
							if (dIntensity > pdTmpRawData[iBinIon]) {
								pdTmpRawData[iBinIon] = dIntensity;
							}
							if (pdTmpRawData[iBinIon] > pPre->dHighestIntensity) {
								pPre->dHighestIntensity = pdTmpRawData[iBinIon];
							}
						}
					} else if (ProNovoConfig::options.iRemovePrecursor == 2) {
						int j;
						int bNotPrec = 1;
						for (j = 1; j <= pScoring->_spectrumInfoInternal.iChargeState; j++) {
							double dMZ;
							dMZ = (mstSpectrum->dParentMass + (j - 1) * PROTON_MASS) / (double) (j);
							if (fabs(dIon - dMZ) < ProNovoConfig::options.dRemovePrecursorTol) {
								bNotPrec = 0;
								break;
							}
						}
						if (bNotPrec) {
							if (dIntensity > pdTmpRawData[iBinIon])
								pdTmpRawData[iBinIon] = dIntensity;

							if (pdTmpRawData[iBinIon] > pPre->dHighestIntensity)
								pPre->dHighestIntensity = pdTmpRawData[iBinIon];
						}
					} else // iRemovePrecursor==0
					{
						if (dIntensity > pdTmpRawData[iBinIon])
							pdTmpRawData[iBinIon] = dIntensity;

						if (pdTmpRawData[iBinIon] > pPre->dHighestIntensity)
							pPre->dHighestIntensity = pdTmpRawData[iBinIon];
					}
				}
			}
		}
	}
	return true;
}

// pdTmpRawData now holds raw data, pdTmpCorrelationData is windowed data after this function
void CometSearchMod::MakeCorrData(double *pdTmpRawData, double *pdTmpCorrelationData, struct Query *pScoring, struct PreprocessStruct *pPre) {
	int i, iWindowSize, iNumWindows = 10;
	double dMaxWindowInten, dTmp1, dTmp2;

	iWindowSize = (int) ((pPre->iHighestIon) / iNumWindows) + 1;
	const int arraySize = pScoring->_spectrumInfoInternal.iArraySize;

	for (i = 0; i < iNumWindows; i++) {
		dMaxWindowInten = 0.0;
		const int windowBegin = i * iWindowSize;
		const int windowEnd = std::min(windowBegin + iWindowSize, arraySize);

#ifndef _MSC_VER
#pragma omp simd reduction(max : dMaxWindowInten)
#endif
		for (int iBin = windowBegin; iBin < windowEnd; ++iBin) {
			if (pdTmpRawData[iBin] > dMaxWindowInten) {
				dMaxWindowInten = pdTmpRawData[iBin];
			}
		}

		if (dMaxWindowInten > 0.0) {
			dTmp1 = 50.0 / dMaxWindowInten;
			dTmp2 = 0.05 * pPre->dHighestIntensity;


#ifndef _MSC_VER
#pragma omp simd
#endif
			for (int iBin = windowBegin; iBin < windowEnd; ++iBin) {
				if (pdTmpRawData[iBin] > dTmp2) {
					pdTmpCorrelationData[iBin] = pdTmpRawData[iBin] * dTmp1;
				}
			}
		}
	}
}

// Smooth input data over 5 points.
bool CometSearchMod::Smooth(double *data, int iArraySize, double *pdTmpSmoothedSpectrum) {
	int i;

	data[0] = 0.0;
	data[1] = 0.0;
	data[iArraySize - 1] = 0.0;
	data[iArraySize - 2] = 0.0;

	for (i = 2; i < iArraySize - 2; i++) {
		// *0.0625 is same as divide by 16.
		pdTmpSmoothedSpectrum[i] = (data[i - 2] + 4.0 * data[i - 1] + 6.0 * data[i] + 4.0 * data[i + 1] + data[i + 2]) * 0.0625;
	}

	memcpy(data, pdTmpSmoothedSpectrum, iArraySize * sizeof(double));

	return true;
}

// Run 2 passes through to pull out peaks.
bool CometSearchMod::PeakExtract(double *data, int iArraySize, double *pdTmpPeakExtracted) {
	int i, ii, iStartIndex, iEndIndex;
	double dStdDev, dAvgInten;

	// 1st pass, choose only peak greater than avg + dStdDev.
	for (i = 0; i < iArraySize; i++) {
		pdTmpPeakExtracted[i] = 0.0;
		dAvgInten = 0.0;

		iStartIndex = i - 50;
		if (i - 50 < 0)
			iStartIndex = 0;

		iEndIndex = i + 50;
		if (i + 50 > iArraySize - 1)
			iEndIndex = iArraySize - 1;

		for (ii = iStartIndex; ii <= iEndIndex; ii++)
			dAvgInten += (double) data[ii];
		dAvgInten /= iEndIndex - iStartIndex;

		dStdDev = 0.0;
		for (ii = iStartIndex; ii <= iEndIndex; ii++)
			dStdDev += (data[ii] - dAvgInten) * (data[ii] - dAvgInten);
		dStdDev = sqrt(dStdDev / (iEndIndex - iStartIndex + 1));

		if ((i > 0) && (i < iArraySize - 1)) {
			if (data[i] > (dAvgInten + dStdDev)) {
				pdTmpPeakExtracted[i] = data[i] - dAvgInten + dStdDev;
				data[i] = 0;     // Remove the peak before 2nd pass.
			}
		}
	}

	// 2nd pass, choose only peak greater than avg + 2*dStdDev.
	for (i = 0; i < iArraySize; i++) {
		dAvgInten = 0.0;

		iStartIndex = i - 50;
		if (i - 50 < 0)
			iStartIndex = 0;

		iEndIndex = i + 50;
		if (i + 50 > iArraySize - 1)
			iEndIndex = iArraySize - 1;

		for (ii = iStartIndex; ii <= iEndIndex; ii++)
			dAvgInten += (double) data[ii];
		dAvgInten /= iEndIndex - iStartIndex;

		dStdDev = 0.0;
		for (ii = iStartIndex; ii <= iEndIndex; ii++)
			dStdDev += (data[ii] - dAvgInten) * (data[ii] - dAvgInten);
		dStdDev = sqrt(dStdDev / (iEndIndex - iStartIndex + 1));

		if ((i > 0) && (i < iArraySize - 1)) {
			if (data[i] > (dAvgInten + 2 * dStdDev))
				pdTmpPeakExtracted[i] = data[i] - dAvgInten + dStdDev;
		}
	}

	memcpy(data, pdTmpPeakExtracted, (size_t) iArraySize * sizeof(double));

	return true;
}

// Pull out top # ions for intensity matching in search.
void CometSearchMod::GetTopIons(double *pdTmpRawData, struct msdata *pTmpSpData, int iArraySize) {
	int i, ii, iLowestIntenIndex = 0;
	double dLowestInten = 0.0, dMaxInten = 0.0;

	for (i = 0; i < iArraySize; i++) {
		if (pdTmpRawData[i] > dLowestInten) {
			//seg debug b
			if (iLowestIntenIndex >= NUM_SP_IONS) {
				cout << "Error 1" << endl;
			}
			//seg debug e
			(pTmpSpData + iLowestIntenIndex)->dIntensity = (double) pdTmpRawData[i];
			(pTmpSpData + iLowestIntenIndex)->dIon = (double) i;

			if ((pTmpSpData + iLowestIntenIndex)->dIntensity > dMaxInten)
				dMaxInten = (pTmpSpData + iLowestIntenIndex)->dIntensity;

			dLowestInten = (pTmpSpData + 0)->dIntensity;
			iLowestIntenIndex = 0;

			for (ii = 1; ii < NUM_SP_IONS; ii++) {
				if ((pTmpSpData + ii)->dIntensity < dLowestInten) {
					dLowestInten = (pTmpSpData + ii)->dIntensity;
					iLowestIntenIndex = ii;
				}
			}
		}
	}

	if (dMaxInten > FLOAT_ZERO) {
		for (i = 0; i < NUM_SP_IONS; i++)
			(pTmpSpData + i)->dIntensity = (((pTmpSpData + i)->dIntensity) / dMaxInten) * 100.0;
	}
}

int CometSearchMod::QsortByIon(const void *p0, const void *p1) {
	if (((struct msdata *) p1)->dIon < ((struct msdata *) p0)->dIon)
		return (1);
	else if (((struct msdata *) p1)->dIon > ((struct msdata *) p0)->dIon)
		return (-1);
	else
		return (0);
}

// Works on Sp data.
void CometSearchMod::StairStep(struct msdata *pTmpSpData, double dFragmentBinSize) {
	int i, ii, iii;
	double dMaxInten, dGap;

	i = 0;
	while (i < NUM_SP_IONS - 1) {
		ii = i;
		dMaxInten = (pTmpSpData + i)->dIntensity;
		dGap = 0.0;

		while (dGap <= dFragmentBinSize && ii < NUM_SP_IONS - 1) {
			ii++;
			dGap = (pTmpSpData + ii)->dIon - (pTmpSpData + ii - 1)->dIon;

			// Finds the max intensity for adjacent points.
			if (dGap <= dFragmentBinSize) {
				if ((pTmpSpData + ii)->dIntensity > dMaxInten)
					dMaxInten = (pTmpSpData + ii)->dIntensity;
			}
		}

		// Sets the adjacent points to the dMaxInten.
		for (iii = i; iii < ii; iii++)
			(pTmpSpData + iii)->dIntensity = dMaxInten;

		i = ii;
	}
}

void CometSearchMod::print(Query * pScoring) {
	int i = 0, x = 0, y = 0;
	cout << "\n" << endl;
	for (i = 1; i < pScoring->_spectrumInfoInternal.iArraySize; i++) {
		x = i / SPARSE_MATRIX_SIZE;
		if (pScoring->ppfSparseFastXcorrData[x] == NULL) {
			cout << "0\t";
		} else {
			y = i - (x * SPARSE_MATRIX_SIZE);
			cout << pScoring->ppfSparseFastXcorrData[x][y] << "\t";
		}
	}
	cout << endl;
	for (i = 1; i < pScoring->_spectrumInfoInternal.iArraySize; i++) {
		x = i / SPARSE_MATRIX_SIZE;
		if (pScoring->ppfSparseFastXcorrDataNL[x] == NULL) {
			cout << "0\t";
		} else {
			y = i - (x * SPARSE_MATRIX_SIZE);
			cout << pScoring->ppfSparseFastXcorrDataNL[x][y] << "\t";
		}
	}
	cout << endl;
	for (i = 1; i < pScoring->_spectrumInfoInternal.iArraySize; i++) {
		x = i / SPARSE_MATRIX_SIZE;
		if (pScoring->ppfSparseSpScoreData[x] == NULL) {
			cout << "0\t";
		} else {
			y = i - (x * SPARSE_MATRIX_SIZE);
			cout << pScoring->ppfSparseSpScoreData[x][y] << "\t";
		}
	}
	cout << "\n" << endl;
}

bool CometSearchMod::ScorePeptides(string * currentPeptide, bool *pbDuplFragment, double * _pdAAforward, double * _pdAAreverse, MS2Scan * mstSpectrum,
		int precursorCharge, unsigned int *** _uiBinnedIonMasses, double & dXcorr) {
	double dInverseBinWidth = 0, dOneMinusBinOffset = 0;
	if (mstSpectrum->isMS2HighRes) {
		dInverseBinWidth = ProNovoConfig::dHighResInverseBinWidth;
		dOneMinusBinOffset = ProNovoConfig::dHighResOneMinusBinOffset;
	} else {
		dInverseBinWidth = ProNovoConfig::dLowResInverseBinWidth;
		dOneMinusBinOffset = ProNovoConfig::dLowResOneMinusBinOffset;
	}

	string * sSequence = currentPeptide;
	size_t i = 0;
	size_t k = sSequence->length() - 1;
	char currentPTM = 0;
	size_t iPeptideLength = 0;
	size_t iPos = 0;
	// map<char, double>::iterator iterResidueMonoMass;
	double iterResidueMonoMass;
	for (i = 0; i <= k; ++i) {
		if (isalpha(sSequence->at(i))) {
			iPeptideLength = iPeptideLength + 1;
		}
	}
	const int iLenMinus1 = static_cast<int>(iPeptideLength) - 1;
	//seg debug b
	if (iLenMinus1 >= iMAX_PEPTIDE_LEN) {
		cout << iLenMinus1 << endl;
		cout << "Error 21" << endl;
	}
	//seg debug e
	if ((int) iPeptideLength < ProNovoConfig::getMinPeptideLength()) {
		cerr << "ERROR: Peptide sequence is too short " << sSequence << endl;
		return false;
	}
	size_t j = 0;
	if (sSequence->at(j) != '[') {
		cerr << "ERROR: First character in a peptide sequence must be [." << endl;
		return false;
	}
	double dBion = ProNovoConfig::precalcMasses.dNtermProton;
	double dYion = ProNovoConfig::precalcMasses.dCtermOH2Proton;
	j++;
	if (!isalpha(sSequence->at(j))) {
		currentPTM = sSequence->at(j);
		iterResidueMonoMass = ProNovoConfig::pdAAMassFragment.find(currentPTM);
		if (iterResidueMonoMass == ProNovoConfig::pdAAMassFragment.end()) {
			cerr << "ERROR: cannot find PTM in the built-in chemistry: " << currentPTM << endl;
			exit(1);
			return false;
		}
		dBion += iterResidueMonoMass;
		j++;
	}
	if (sSequence->at(k) == ']') {
		k--;
	} else {
		currentPTM = sSequence->at(k);
		iterResidueMonoMass = ProNovoConfig::pdAAMassFragment.find(currentPTM);
		if (iterResidueMonoMass == ProNovoConfig::pdAAMassFragment.end()) {
			cerr << "ERROR: cannot find PTM in the built-in chemistry: " << currentPTM << endl;
			exit(1);
			return false;
		}
		dYion += iterResidueMonoMass;
		k--;
		if (sSequence->at(k) != ']') {
			cerr << "ERROR: (second) Last character in a peptide sequence must be ]." << endl;
			exit(1);
			return false;
		}
		k--;
	}
	for (i = 0; i < (size_t) iLenMinus1; i++) {
		//First forward
		if (!isalpha(sSequence->at(j))) {
			cerr << "ERROR: One residue can only have one PTM (Up to only one symbol after an amino acid)" << endl;
			exit(1);
			return false;
		}
		currentPTM = sSequence->at(j);
		iterResidueMonoMass = ProNovoConfig::pdAAMassFragment.find(currentPTM);
		if (iterResidueMonoMass == ProNovoConfig::pdAAMassFragment.end()) {
			cerr << "ERROR: One residue can only have one PTM (Up to only one symbol after an amino acid)" << endl;
			exit(1);
			return false;
		}
		dBion += iterResidueMonoMass;
		j++;
		if (!isalpha(sSequence->at(j))) {
			currentPTM = sSequence->at(j);
			iterResidueMonoMass = ProNovoConfig::pdAAMassFragment.find(currentPTM);
			if (iterResidueMonoMass == ProNovoConfig::pdAAMassFragment.end()) {
				cerr << "ERROR: cannot find PTM in the built-in chemistry: " << currentPTM << endl;
				exit(1);
				return false;
			}
			dBion += iterResidueMonoMass;
			j++;
		}
		//seg debug b
		if ((int) iPos >= iMAX_PEPTIDE_LEN) {
			cout << iPos << endl;
			cout << "Error 20" << endl;
		}
		//seg debug e
		_pdAAforward[iPos] = dBion;
		//Now reverse
		if (!isalpha(sSequence->at(k))) {
			currentPTM = sSequence->at(k);
			iterResidueMonoMass = ProNovoConfig::pdAAMassFragment.find(currentPTM);
			if (iterResidueMonoMass == ProNovoConfig::pdAAMassFragment.end()) {
				cerr << "ERROR: cannot find PTM in the built-in chemistry: " << currentPTM << endl;
				exit(1);
				return false;
			}
			dYion += iterResidueMonoMass;
			k--;
		}
		if (!isalpha(sSequence->at(k))) {
			cerr << "ERROR: One residue can only have one PTM (Up to only one symbol after an amino acid)" << endl;
			exit(1);
			return false;
		}
		currentPTM = sSequence->at(k);
		iterResidueMonoMass = ProNovoConfig::pdAAMassFragment.find(currentPTM);
		if (iterResidueMonoMass == ProNovoConfig::pdAAMassFragment.end()) {
			cerr << "ERROR: cannot find PTM in the built-in chemistry: " << currentPTM << endl;
			exit(1);
			return false;
		}
		dYion += iterResidueMonoMass;
		k--;
		//seg debug b
		if ((int) iPos >= iMAX_PEPTIDE_LEN) {
			cout << iPos << endl;
			cout << "Error 20" << endl;
		}
		//seg debug e
		_pdAAreverse[iPos] = dYion;
		iPos++;
	}
	// seg debug b
	bool bxx = false;
	for (int xx = 0; xx < iLenMinus1; xx++) {
		if (_pdAAforward[xx] >= (mstSpectrum->dParentMass + 5.0)) {
			cout << xx << endl;
			cout << _pdAAforward[xx] << endl;
			cout << "Error 50" << endl;
		}
		if (_pdAAreverse[xx] >= (mstSpectrum->dParentMass + 5.0)) {
			cout << xx << endl;
			cout << _pdAAreverse[xx] << endl;
			cout << "Error 51" << endl;
			bxx = true;
		}
	}
	if (bxx) {
		cout << (*sSequence) << endl;
		cout << ProNovoConfig::pdAAMassFragment.find('o') << endl;
		cout << ProNovoConfig::pdAAMassFragment.find('h') << endl;
		dYion = ProNovoConfig::precalcMasses.dCtermOH2Proton;
		cout << dYion << endl;
		k = sSequence->length() - 1;
		if (sSequence->at(k) == ']') {
			k--;
		} else {
			currentPTM = sSequence->at(k);
			cout << k << ": " << currentPTM << endl;
			iterResidueMonoMass = ProNovoConfig::pdAAMassFragment.find(currentPTM);
			if (iterResidueMonoMass == ProNovoConfig::pdAAMassFragment.end()) {
				cerr << "ERROR: cannot find PTM in the built-in chemistry: " << currentPTM << endl;
				return false;
			}
			dYion += iterResidueMonoMass;
			cout << ".. " << iterResidueMonoMass << endl;
			cout << ", " << dYion << endl;
			k--;
			if (sSequence->at(k) != ']') {
				cerr << "ERROR: (second) Last character in a peptide sequence must be ]." << endl;
				return false;
			}
			k--;
		}
		for (i = 0; i < (size_t) iLenMinus1; i++) {
			//Now reverse
			if (!isalpha(sSequence->at(k))) {
				currentPTM = sSequence->at(k);
				cout << k << ": " << currentPTM << endl;
				iterResidueMonoMass = ProNovoConfig::pdAAMassFragment.find(currentPTM);
				if (iterResidueMonoMass == ProNovoConfig::pdAAMassFragment.end()) {
					cerr << "ERROR: cannot find PTM in the built-in chemistry: " << currentPTM << endl;
					return false;
				}
				dYion += iterResidueMonoMass;
				cout << ".. " << iterResidueMonoMass << endl;
				k--;
			}
			if (!isalpha(sSequence->at(k))) {
				cerr << "ERROR: One residue can only have one PTM (Up to only one symbol after an amino acid)" << endl;
				return false;
			}
			currentPTM = sSequence->at(k);
			cout << k << ": " << currentPTM << endl;
			iterResidueMonoMass = ProNovoConfig::pdAAMassFragment.find(currentPTM);
			if (iterResidueMonoMass == ProNovoConfig::pdAAMassFragment.end()) {
				cerr << "ERROR: cannot find PTM in the built-in chemistry: " << currentPTM << endl;
				return false;
			}
			dYion += iterResidueMonoMass;
			cout << ".. " << iterResidueMonoMass << endl;
			k--;
			cout << i << ", " << dYion << endl;
		}
	}
	// seg debug e
	// Now get the set of binned fragment ions once to compare this peptide against all matching spectra.
	int ctCharge = 0;
	int ctIonSeries = 0;
	int iWhichIonSeries = 0;
	int ctLen = 0;
	int iVal = 0;
	Query * pQuery = mstSpectrum->pQuery;
	// seg debug b
	if (pQuery == NULL) {
		cout << "Error 96" << endl;
	}
	// seg debug e
	int iMaxFragCharge = pQuery->_spectrumInfoInternal.iMaxFragCharge;
	if (!mstSpectrum->bUseReactionChargeForScoring) {
		iMaxFragCharge = precursorCharge <= 1 ? 1 : precursorCharge - 1;
	}
	// seg debug b
	if (iMaxFragCharge >= iMaxPercusorCharge) {
		cout << "Error 9" << endl;
	}
	if (pbDuplFragment == NULL) {
		cout << "Error 98" << endl;
	}
	// seg debug e
	for (ctCharge = 1; ctCharge <= iMaxFragCharge; ctCharge++) {
		for (ctIonSeries = 0; ctIonSeries < ProNovoConfig::ionInformation.iNumIonSeriesUsed; ctIonSeries++) {
			iWhichIonSeries = ProNovoConfig::ionInformation.piSelectedIonSeries[ctIonSeries];
			// seg debug b
			if (iWhichIonSeries >= iDimesion2) {
				cout << "Error 10" << endl;
			}
			// seg debug e
			for (ctLen = 0; ctLen < iLenMinus1; ctLen++) {
				iVal = BINX(GetFragmentIonMass(iWhichIonSeries, ctLen, ctCharge, _pdAAforward, _pdAAreverse), dInverseBinWidth, dOneMinusBinOffset);
				//seg debug b
				if (iVal >= iArraySizeScore) {
					cout << iArraySizeScore << endl;
					cout << iVal << endl;
					cout << ProNovoConfig::dMaxPeptideMass << endl;
					cout << "Error 7" << endl;
				}
				if (iVal >= pQuery->_spectrumInfoInternal.iArraySize) {
					cout << pQuery->_spectrumInfoInternal.iArraySize << endl;
					cout << iVal << endl;
					cout << ctLen << endl;
					cout << iWhichIonSeries << endl;
					cout << _pdAAforward[ctLen] << endl;
					cout << _pdAAreverse[ctLen] << endl;
					cout << mstSpectrum->dParentMass << endl;
					cout << "Error 30" << endl;
					exit(1);
				}
				//seg debug e
				pbDuplFragment[iVal] = false;
			}
		}
	}
	// seg debug b
	if (_uiBinnedIonMasses == NULL) {
		cout << "Error 97" << endl;
	}
	// seg debug e
	for (ctCharge = 1; ctCharge <= iMaxFragCharge; ctCharge++) {
		for (ctIonSeries = 0; ctIonSeries < ProNovoConfig::ionInformation.iNumIonSeriesUsed; ctIonSeries++) {
			iWhichIonSeries = ProNovoConfig::ionInformation.piSelectedIonSeries[ctIonSeries];
			// seg debug b
			if (iWhichIonSeries >= iDimesion2) {
				cout << "Error 10" << endl;
			}
			// seg debug e
			for (ctLen = 0; ctLen < iLenMinus1; ctLen++) {
				iVal = BINX(GetFragmentIonMass(iWhichIonSeries, ctLen, ctCharge, _pdAAforward, _pdAAreverse), dInverseBinWidth, dOneMinusBinOffset);
				//seg debug b
				if (iVal >= iArraySizeScore) {
					cout << iArraySizeScore << endl;
					cout << iVal << endl;
					cout << "Error 8" << endl;
				}
				if (iVal >= pQuery->_spectrumInfoInternal.iArraySize) {
					cout << pQuery->_spectrumInfoInternal.iArraySize << endl;
					cout << iVal << endl;
					cout << "Error 31" << endl;
				}
				//seg debug e
				if (pbDuplFragment[iVal] == false) {
					_uiBinnedIonMasses[ctCharge][ctIonSeries][ctLen] = iVal;
					pbDuplFragment[iVal] = true;
				} else {
					_uiBinnedIonMasses[ctCharge][ctIonSeries][ctLen] = 0;
				}
			}
		}
	}

	dXcorr = 0;
	bool bUseNLPeaks = false;
	float **ppSparseFastXcorrData;              // use this if bSparseMatrix
	int bin, x, y;
	int iMax = pQuery->_spectrumInfoInternal.iArraySize / SPARSE_MATRIX_SIZE + 1;
	for (ctCharge = 1; ctCharge <= iMaxFragCharge; ctCharge++) {
		for (ctIonSeries = 0; ctIonSeries < ProNovoConfig::ionInformation.iNumIonSeriesUsed; ctIonSeries++) {
			iWhichIonSeries = ProNovoConfig::ionInformation.piSelectedIonSeries[ctIonSeries];
			if (ProNovoConfig::ionInformation.bUseNeutralLoss
					&& (ProNovoConfig::ionInformation.iIonVal[ION_SERIES_A] || ProNovoConfig::ionInformation.iIonVal[ION_SERIES_B]
							|| ProNovoConfig::ionInformation.iIonVal[ION_SERIES_Y])) {
				bUseNLPeaks = true;
			} else {
				bUseNLPeaks = false;
			}
			if (ctCharge == 1 && bUseNLPeaks) {
				ppSparseFastXcorrData = pQuery->ppfSparseFastXcorrDataNL;
			} else {
				ppSparseFastXcorrData = pQuery->ppfSparseFastXcorrData;
			}
			// seg debug b
			if (ppSparseFastXcorrData == NULL) {
				cout << "Error 98" << endl;
			}
			// seg debug e
			for (ctLen = 0; ctLen < iLenMinus1; ctLen++) {
				//MH: newer sparse matrix converts bin to sparse matrix bin
				bin = _uiBinnedIonMasses[ctCharge][ctIonSeries][ctLen];
				// seg debug b
				if (bin >= pQuery->_spectrumInfoInternal.iArraySize) {
					cout << "Error 32" << endl;
				}
				// seg debug e
				x = bin / SPARSE_MATRIX_SIZE;
				if (ppSparseFastXcorrData[x] == NULL || x > iMax) // x should never be > iMax so this is just a safety check
					continue;
				y = bin - (x * SPARSE_MATRIX_SIZE);
				dXcorr += ppSparseFastXcorrData[x][y];
				// cout << bin << endl;
				//cout << x << "," << y << endl;
			}
		}
	}
	if (dXcorr < XCORR_CUTOFF) {
		dXcorr = XCORR_CUTOFF;
		return false;
	} else {
		dXcorr *= 0.005;  // Scale intensities to 50 and divide score by 1E4.
	}

	/*cout << "\n" << endl;
	 for (i = 1; (int) i < pQuery->_spectrumInfoInternal.iArraySize; i++) {
	 x = i / SPARSE_MATRIX_SIZE;
	 if (pQuery->ppfSparseFastXcorrDataNL[x] == NULL) {
	 //cout << "0\t";
	 } else {
	 y = i - (x * SPARSE_MATRIX_SIZE);
	 if (pQuery->ppfSparseFastXcorrDataNL[x][y] != 0) {
	 cout << pQuery->ppfSparseFastXcorrDataNL[x][y] << "\t";
	 }
	 }
	 }
	 cout << endl;
	 for (i = 1; (int) i < pQuery->_spectrumInfoInternal.iArraySize; i++) {
	 x = i / SPARSE_MATRIX_SIZE;
	 if (pQuery->ppfSparseFastXcorrData[x] == NULL) {
	 //cout << "0\t";
	 } else {
	 y = i - (x * SPARSE_MATRIX_SIZE);
	 if (pQuery->ppfSparseFastXcorrData[x][y] != 0) {
	 cout << pQuery->ppfSparseFastXcorrData[x][y] << "\t";
	 }
	 }
	 }
	 cout << endl;*/

	return true;
}

void CometSearchMod::BinPeptideIonsSIPNoCancelOut(
		const vector<vector<double> > &vvdYionMass,
		const vector<vector<double> > &vvdYionProb,
		const vector<vector<double> > &vvdBionMass,
		const vector<vector<double> > &vvdBionProb,
		MS2Scan *mstSpectrum,
		vector<unsigned char> &pbDuplFragment,
		vector<double> &vdBinnedIonMasses,
		vector<int> &vdBin) {
	double dInverseBinWidth = 0, dOneMinusBinOffset = 0;
	// Reuse the sparse list from the preceding candidate to reset only bins
	// that were touched.  The old implementation visited every theoretical
	// ion once to clear its marker, visited all ions again to populate weights,
	// and then allocated a std::set to deduplicate bins.
	for (int bin : vdBin) {
		pbDuplFragment.at(bin) = 0;
		vdBinnedIonMasses.at(bin) = 0.0;
	}
	vdBin.clear();
	if (mstSpectrum->isMS2HighRes) {
		dInverseBinWidth = ProNovoConfig::dHighResInverseBinWidth;
		dOneMinusBinOffset = ProNovoConfig::dHighResOneMinusBinOffset;
	} else {
		dInverseBinWidth = ProNovoConfig::dLowResInverseBinWidth;
		dOneMinusBinOffset = ProNovoConfig::dLowResOneMinusBinOffset;
	}

	// Now get the set of binned fragment ions once to compare this peptide against all matching spectra.
	int ctCharge = 0;
	// int ctIonSeries = 0;
	// int iWhichIonSeries = 0;
	// int ctLen = 0;
	int iVal = 0;
	Query * pQuery = mstSpectrum->pQuery;

	int iMaxFragCharge = pQuery->_spectrumInfoInternal.iMaxFragCharge;

	int iIonMassSize = 0, iIsotopicDistSize = 0;
	double dFragmentIonMassZ = 0;
	for (ctCharge = 1; ctCharge <= iMaxFragCharge; ctCharge++) {
		// y-ion
		iIonMassSize = static_cast<int>(vvdYionMass.size());
		for (int i = 0; i < iIonMassSize; ++i) {
			iIsotopicDistSize = static_cast<int>(vvdYionMass.at(i).size());
			for (int j = 0; j < iIsotopicDistSize; ++j) {
				if (vvdYionProb.at(i).at(j) < ProbabilityCutOff) {
					continue;
				}
				dFragmentIonMassZ = (vvdYionMass.at(i).at(j) + (ctCharge) * PROTON_MASS) / ctCharge;
				iVal = BINX(dFragmentIonMassZ, dInverseBinWidth, dOneMinusBinOffset);
				if (iVal < 0 ||
					static_cast<size_t>(iVal) >= vdBinnedIonMasses.size()) {
					continue;
				}
				if (pbDuplFragment.at(iVal) == 0) {
					vdBinnedIonMasses.at(iVal) = vvdYionProb.at(i).at(j);
					pbDuplFragment.at(iVal) = 1;
					vdBin.push_back(iVal);
				} else {
					if (vdBinnedIonMasses.at(iVal) < vvdYionProb.at(i).at(j)) {
						vdBinnedIonMasses.at(iVal) = vvdYionProb.at(i).at(j);
					}
				}
			}
		}
		// b-ion
		iIonMassSize = static_cast<int>(vvdBionMass.size());
		for (int i = 0; i < iIonMassSize; ++i) {
			iIsotopicDistSize = static_cast<int>(vvdBionMass.at(i).size());
			for (int j = 0; j < iIsotopicDistSize; ++j) {
				if (vvdBionProb.at(i).at(j) < ProbabilityCutOff) {
					continue;
				}
				dFragmentIonMassZ = (vvdBionMass.at(i).at(j) + (ctCharge) * PROTON_MASS) / ctCharge;
				iVal = BINX(dFragmentIonMassZ, dInverseBinWidth, dOneMinusBinOffset);
				if (iVal < 0 ||
					static_cast<size_t>(iVal) >= vdBinnedIonMasses.size()) {
					continue;
				}
				if (pbDuplFragment.at(iVal) == 0) {
					vdBinnedIonMasses.at(iVal) = vvdBionProb.at(i).at(j);
					pbDuplFragment.at(iVal) = 1;
					vdBin.push_back(iVal);
				} else {
					if (vdBinnedIonMasses.at(iVal) < vvdBionProb.at(i).at(j)) {
						vdBinnedIonMasses.at(iVal) = vvdBionProb.at(i).at(j);
					}
				}
			}
		}
	}

	std::sort(vdBin.begin(), vdBin.end());

	vdBinnedIonMasses.at(0) = 1;
}

bool CometSearchMod::ScoreBinnedPeptideSIPNoCancelOut(
		MS2Scan *mstSpectrum,
		const vector<double> &vdBinnedIonMasses,
		const vector<int> &vdBin,
		double &dXcorr) {
	Query *pQuery = mstSpectrum->pQuery;
	dXcorr = 0;
	// bool bUseNLPeaks = false;
	float **ppSparseFastXcorrData;              // use this if bSparseMatrix
	int bin, x, y;
	int iMax = pQuery->_spectrumInfoInternal.iArraySize / SPARSE_MATRIX_SIZE + 1;
	ppSparseFastXcorrData = pQuery->ppfSparseFastXcorrData;
	if (ppSparseFastXcorrData == NULL) {
		cout << "error 20170308" << endl;
		dXcorr = XCORR_CUTOFF;
		return false;
	}
	const int iBinSize = static_cast<int>(vdBin.size());
	for (int i = 0; i < iBinSize; ++i) {
		bin = vdBin.at(i);
		if (bin < 0 || bin >= pQuery->_spectrumInfoInternal.iArraySize)
			continue;
		x = bin / SPARSE_MATRIX_SIZE;
		if (x > iMax || ppSparseFastXcorrData[x] == NULL)
			continue;
		y = bin - (x * SPARSE_MATRIX_SIZE);
		dXcorr += ppSparseFastXcorrData[x][y] * vdBinnedIonMasses.at(bin);
		// cout << bin << endl;
	}
	if (dXcorr < XCORR_CUTOFF) {
		dXcorr = XCORR_CUTOFF;
		return false;
	} else {
		dXcorr *= 0.005;  // Scale intensities to 50 and divide score by 1E4.
	}

	return true;
}

bool CometSearchMod::ScorePeptidesSIPNoCancelOut(const vector<vector<double> > & vvdYionMass, const vector<vector<double> > & vvdYionProb,
		const vector<vector<double> > & vvdBionMass, const vector<vector<double> > & vvdBionProb, MS2Scan * mstSpectrum, vector<unsigned char> & pbDuplFragment,
		vector<double> & vdBinnedIonMasses, vector<int> & vdBin, double & dXcorr) {
	BinPeptideIonsSIPNoCancelOut(
		vvdYionMass, vvdYionProb, vvdBionMass, vvdBionProb,
		mstSpectrum, pbDuplFragment, vdBinnedIonMasses, vdBin);
	return ScoreBinnedPeptideSIPNoCancelOut(
		mstSpectrum, vdBinnedIonMasses, vdBin, dXcorr);
}

double CometSearchMod::GetFragmentIonMass(int iWhichIonSeries, int i, int ctCharge, double *_pdAAforward, double *_pdAAreverse) {
	double dFragmentIonMass = 0.0;

	switch (iWhichIonSeries) {
	case ION_SERIES_B:
		dFragmentIonMass = _pdAAforward[i];
		break;
	case ION_SERIES_Y:
		dFragmentIonMass = _pdAAreverse[i];
		break;
	case ION_SERIES_A:
		cout << "Error 40" << endl;
		dFragmentIonMass = _pdAAforward[i] - ProNovoConfig::precalcMasses.dCO;
		break;
	case ION_SERIES_C:
		cout << "Error 41" << endl;
		dFragmentIonMass = _pdAAforward[i] + ProNovoConfig::precalcMasses.dNH3;
		break;
	case ION_SERIES_Z:
		cout << "Error 42" << endl;
		dFragmentIonMass = _pdAAreverse[i] - ProNovoConfig::precalcMasses.dNH2;
		break;
	case ION_SERIES_X:
		cout << "Error 43" << endl;
		dFragmentIonMass = _pdAAreverse[i] + ProNovoConfig::precalcMasses.dCOminusH2;
		break;
	}

	return (dFragmentIonMass + (ctCharge - 1) * PROTON_MASS) / ctCharge;
}

bool CometSearchMod::CalculateSP(double & fScoreSp, double* _pdAAforward, double * _pdAAreverse, MS2Scan * mstSpectrum, int iLenPeptide) {
	double dInverseBinWidth = 0, dOneMinusBinOffset = 0;
	if (mstSpectrum->isMS2HighRes) {
		dInverseBinWidth = ProNovoConfig::dHighResInverseBinWidth;
		dOneMinusBinOffset = ProNovoConfig::dHighResOneMinusBinOffset;
	} else {
		dInverseBinWidth = ProNovoConfig::dLowResInverseBinWidth;
		dOneMinusBinOffset = ProNovoConfig::dLowResOneMinusBinOffset;
	}
	int iMaxFragCharge = mstSpectrum->pQuery->_spectrumInfoInternal.iMaxFragCharge;
	int iMatchedFragmentIonCt = 0;
	IonSeriesStruct ionSeries[9];
	double dTmpIntenMatch = 0.0;
	double dConsec = 0.0;
	int ii;
	for (ii = 0; ii < ProNovoConfig::ionInformation.iNumIonSeriesUsed; ii++) {
		int iii;
		for (iii = 1; iii <= iMaxFragCharge; iii++)
			ionSeries[ProNovoConfig::ionInformation.piSelectedIonSeries[ii]].bPreviousMatch[iii] = 0;
	}
	int ctCharge;
	double dFragmentIonMass = 0.0;
	for (ctCharge = 1; ctCharge <= iMaxFragCharge; ctCharge++) {
		for (ii = 0; ii < ProNovoConfig::ionInformation.iNumIonSeriesUsed; ii++) {
			int iWhichIonSeries = ProNovoConfig::ionInformation.piSelectedIonSeries[ii];

			// As both _pdAAforward and _pdAAreverse are increasing, loop through
			// iLenPeptide-1 to complete set of internal fragment ions.
			for (int iii = 0; iii < iLenPeptide - 1; iii++) {
				// Gets fragment ion mass.
				dFragmentIonMass = GetFragmentIonMass(iWhichIonSeries, iii, ctCharge, _pdAAforward, _pdAAreverse);

				if (!(dFragmentIonMass <= FLOAT_ZERO)) {
					int iFragmentIonMass = BINX(dFragmentIonMass, dInverseBinWidth, dOneMinusBinOffset);
					float fSpScore;

					fSpScore = static_cast<float>(
						FindSpScore(mstSpectrum->pQuery, iFragmentIonMass));

					if (fSpScore > FLOAT_ZERO) {
						iMatchedFragmentIonCt++;

						// Simple sum intensity.
						dTmpIntenMatch += fSpScore;

						// Increase score for consecutive fragment ion series matches.
						if (ionSeries[iWhichIonSeries].bPreviousMatch[ctCharge])
							dConsec += 0.075;

						ionSeries[iWhichIonSeries].bPreviousMatch[ctCharge] = 1;
					} else {
						ionSeries[iWhichIonSeries].bPreviousMatch[ctCharge] = 0;
					}
				}
			}
		}
	}
	fScoreSp = (double) ((dTmpIntenMatch * iMatchedFragmentIonCt * (1.0 + dConsec))
			/ ((iLenPeptide) * iMaxFragCharge * ProNovoConfig::ionInformation.iNumIonSeriesUsed));

	return true;
}

double CometSearchMod::FindSpScore(Query *pQuery, int bin) {
	int x = bin / SPARSE_MATRIX_SIZE;
	if (pQuery->ppfSparseSpScoreData[x] == NULL)
		return 0.0f;
	int y = bin - (x * SPARSE_MATRIX_SIZE);
	return pQuery->ppfSparseSpScoreData[x][y];
}
