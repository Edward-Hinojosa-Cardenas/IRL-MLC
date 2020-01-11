import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;
import java.util.ArrayList;
import java.util.Random;

public class CMain {
	public static void main(String[] args) throws FileNotFoundException, IOException {
		String sDS = "scene";
		int iCV = 10;
		int iFold = 10;
		int nAttr = 294;
		int nLabe = 6;
		
		//Parameters of the genetic learning of preliminary rule base
		int nFS = 5;
		int rbs = 400;
		int ps = 200;
		double pdc = 0.9;
		int ni = 500;
		int ts = 2;
		double mp = 0.2;
		double mp_a = 0.2;
		
		//Parameters of the genetic tuning of preliminary rule base
		int pst = 100;
		int nit = 1000;
		int tst = 2;
		double dAlfa = 0.5;
		double mpt = 0.1;
		double mpt_wth = 0.1;
		
		System.out.println("Dataset = " + sDS);
		System.out.println("Cross validation = " + iCV);
		System.out.println("Cross validation fold= " + iFold);
		System.out.println("Number of attributes = " + nAttr);
		System.out.println("Number of labels = " + nLabe);
		
		System.out.println("************* Parameters of the genetic learning of preliminary rule base *******");
		System.out.println("Number of fuzzy sets (nFS) = " + nFS);
		System.out.println("Preliminar rule base size (rbs) = " + rbs);
		System.out.println("Population size (ps) = " + ps);
		System.out.println("Probability of don't care conditions (pdc) = " + pdc);
		System.out.println("Number of iterations of genetic algorithm (ni) = " + ni);
		System.out.println("Number of fuzzy rules in tournament selection (ts) = " + ts);
		System.out.println("Mutation probability (mp) = " + mp);
		
		System.out.println("************* Parameters of the genetic tuning of preliminary rule base *******");
		System.out.println("Population size (pst) = " + pst);
		System.out.println("Number of iterations of genetic algorithm (nit) = " + nit);
		System.out.println("Number of fuzzy rule weights in tournament selection (tst) = " + tst);
		System.out.println("SBX crossover (Alpha) = " + dAlfa);
		System.out.println("Mutation probability (mpt) = " + mpt);
		System.out.println("Mutarion probability for each weight and threrhold (mpt_wth) = " + mpt_wth);
		System.out.println("*************");
		
		System.out.println();

		int nExamTrain = 0;
		int nExamTest = 0;
		double dataTrain[][] = null;
		double dataTest[][] = null;

		int r = 0;
		double dCantClas = 0;
		
		// Genetic Algorithm - Parameters
		int i1, i2;
		int iMax1, iMax2, iMaxi;
		double dAccu1, dAccu2;
		
		// Read Train and Test Data
		String docuTrain = System.getProperty("user.dir") + File.separator + "datasets" + File.separator
				+ sDS + "-" + iCV + "-fold" + File.separator + sDS + "-" + iCV + "-" + iFold + "tra.dat";

		String docuTest = System.getProperty("user.dir") + File.separator + "datasets" + File.separator
				+ sDS + "-" + iCV + "-fold" + File.separator + sDS + "-" + iCV + "-" + iFold + "tst.dat";	

		//Read train dataset
		try (BufferedReader br = Files.newBufferedReader(Paths.get(docuTrain))) {
			boolean bData = false;
			String sLine;
			while ((sLine = br.readLine()) != null) {
				if (bData) {
					nExamTrain++;
				}
				if (sLine.contains("@data")) {
					bData = true;
				}
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
		int aRealClasTrain[][] = new int[nExamTrain][nLabe];
		try (BufferedReader br = Files.newBufferedReader(Paths.get(docuTrain))) {
			String sLine;
			String[] asData;

			boolean bData = false;
			r = 0;
			while ((sLine = br.readLine()) != null) {
				if (bData) {
					dCantClas = 0;
					asData = sLine.split(",");
					for (int c = 0; c < asData.length - nLabe; c++) {
						dataTrain[r][c] = Double.parseDouble(asData[c]);
					}
					for (int c = nAttr; c < asData.length; c++) {
						if(asData[c].equals("1")) {
							dCantClas++;
							aRealClasTrain[r][c - nAttr] = 1;
						}
					}
					for (int c = nAttr; c < asData.length; c++) {
						if(asData[c].equals("1")) {
							dataTrain[r][c] = 1.0 / dCantClas;
						}
					}
					r++;
				}
				if (sLine.contains("@data")) {
					bData = true;
					dataTrain = new double[nExamTrain][nAttr + nLabe];
				}
			}
		} catch (IOException e) {
			System.err.format("IOException: %s%n", e);
		}

		//Read test dataset
		try (BufferedReader br = Files.newBufferedReader(Paths.get(docuTest))) {
			boolean bData = false;
			String sLine;
			while ((sLine = br.readLine()) != null) {
				if (bData) {
					nExamTest++;
				}
				if (sLine.contains("@data")) {
					bData = true;
				}
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
		int aRealClasTest[][] = new int[nExamTest][nLabe];
		try (BufferedReader br = Files.newBufferedReader(Paths.get(docuTest))) {
			String sLine;
			String[] asData;

			boolean bData = false;
			r = 0;
			while ((sLine = br.readLine()) != null) {
				if (bData) {
					asData = sLine.split(",");
					for (int c = 0; c < asData.length - nLabe; c++) {
						dataTest[r][c] = Double.parseDouble(asData[c]);
					}
					for (int c = nAttr; c < asData.length; c++) {
						if(asData[c].equals("1")) {
							aRealClasTest[r][c - nAttr] = 1;
						}
					}
					r++;
				}
				if (sLine.contains("@data")) {
					bData = true;
					dataTest = new double[nExamTest][nAttr + nLabe];
				}
			}
		} catch (IOException e) {
			System.err.format("IOException: %s%n", e);
		}
		//showMatrizDouble(dataTrain);
		////////////////////////////////////////////////////////////////////////////////////

		// Create Uniformly Fuzzy Sets
		double dMax;
		double dMin;
		double dStep;
		double dFuzzySets[][] = new double[nAttr][nFS * 3];
		int r_aux;
		int c_aux;

		for (int c = 0; c < nAttr; c++) {
			dMax = Double.MAX_VALUE * -1;
			dMin = Double.MAX_VALUE;
			for (r = 0; r < dataTrain.length; r++) {
				if (dataTrain[r][c] > dMax) {
					dMax = dataTrain[r][c];
				}
				if (dataTrain[r][c] < dMin) {
					dMin = dataTrain[r][c];
				}
			}
			r_aux = c;
			dStep = (dMax - dMin) / (nFS - 1);
			dMin = dMin - dStep;
			for (c_aux = 0; c_aux < nFS * 3; c_aux = c_aux + 3) {
				dFuzzySets[r_aux][c_aux] = dMin;
				dFuzzySets[r_aux][c_aux + 1] = dMin + dStep;
				dFuzzySets[r_aux][c_aux + 2] = dMin + dStep + dStep;
				dMin = dMin + dStep;
			}
		}

		//showMatrizDouble(dFuzzySets);

		// Create Fuzzy Rules
		int sizeRule = nAttr + 2;
		ArrayList<double[]> popuLear_P_O = new ArrayList<>();
		ArrayList<CChrom> preliminaryRB = new ArrayList<>();

		double dHammTrain = -1;
		double dHammTest = -1;

		int iData;
		int iFZ = -1;
		double dMembDegr;
		double dMinMembDegr;
		ArrayList<Integer> iDataSele = new ArrayList<>();
		int iO = 0;
		
		while(preliminaryRB.size() < rbs) {
			// Create Initial Population
			System.out.println("Genetic learning of fuzzy rule " + (preliminaryRB.size() + 1) + " of " + rbs);
			int j = 0;
			iDataSele.clear();
			double[] rule;
			while(j < ps) {
				rule = new double[sizeRule];
				iData = new Random().nextInt(dataTrain.length);
				for (int posRule = 0; posRule < nAttr; posRule++) {
					if (new Random().nextDouble() > pdc) {
						dMinMembDegr = 0.0;
						iFZ = 0;
						for(int x = 1; x < nFS + 1; x++) {
							dMembDegr = getMembershipGrade(dFuzzySets, dataTrain[iData][posRule], posRule, x);
							if(dMembDegr > dMinMembDegr) {
								dMinMembDegr = dMembDegr;
								iFZ = x;
							}
						}
					}
					else {
						iFZ = 0;
					}
					rule[posRule] = iFZ;
				}
				int iInt = 0;
				while(true) {
					iO = new Random().nextInt(nLabe);
					if (dataTrain[iData][nAttr + iO] > 0.0) {
						rule[sizeRule - 2] = iO;
						break;
					}
					iInt++;
					if(iInt >= 50) {
						break;
					}
				}
				rule[sizeRule - 1] = 0.0;
				popuLear_P_O.add(rule);
				//showRule(rule);
				j++;
			}
			//showRB(trainRB, Iinputs, 1, sizeRule, 0);
			// Calculate the fitness of genetic learning
			double dFitnLear;
			for (int i = 0; i < ps; i++) {
				dFitnLear = getFitnessLearning(popuLear_P_O.get(i), dataTrain, dFuzzySets, nAttr, nLabe, sizeRule);
				popuLear_P_O.get(i)[sizeRule - 1] = dFitnLear;
				//System.out.println(dAccu);
			}
			//Iterations
			for (int iT = 0; iT < ni; iT++) {
				//System.out.println("\tIteration of Genetic Algorithm " + (iT + 1) + " of " + iIter);
				for (int iI = 0; iI < ps; iI = iI + 2) {
					// Tournament Selection
					iMax1 = iMax2 = iMaxi = -1;
					for(int t = 0; t < ts; t++) {
						i1 = new Random().nextInt(ps);
						i2 = new Random().nextInt(ps);
						while (i1 == i2) {
							i2 = new Random().nextInt(ps);
						}
						dAccu1 = popuLear_P_O.get(i1)[sizeRule - 1];
						dAccu2 = popuLear_P_O.get(i2)[sizeRule - 1];
						if (dAccu1 > dAccu2) {
							iMaxi = i1;
						} else if (dAccu2 > dAccu1) {
							iMaxi = i2;
						} else {
							if (new Random().nextDouble() < 0.5) {
								iMaxi = i1;
							} else {
								iMaxi = i2;
							}
						}
						if(t == 0) {
							iMax1 = iMaxi;
						}
						else {
							iMax2 = iMaxi;
						}
					}
					double[] offs1 = new double[sizeRule];
					double[] offs2 = new double[sizeRule];
					// Crossover
					double nRand;
					for (j = 0; j < sizeRule - 1; j++) {
						nRand = new Random().nextDouble();
						if (nRand < 0.5) {
							offs1[j] = popuLear_P_O.get(iMax1)[j];
							offs2[j] = popuLear_P_O.get(iMax2)[j];
						}
						else {
							offs1[j] = popuLear_P_O.get(iMax2)[j];
							offs2[j] = popuLear_P_O.get(iMax1)[j];
						}
					}

					// Mutation Descendant 1
					if(new Random().nextDouble() < mp) {
						for (j = 0; j < sizeRule - 2; j++) {
							if(new Random().nextDouble() < mp_a) {
								offs1[j] = new Random().nextInt(nFS + 1);
							}
						}
					}
					// Mutation Descendant 2
					if(new Random().nextDouble() < mp) {
						for (j = 0; j < sizeRule - 2; j++) {
							if(new Random().nextDouble() < mp_a) {
								offs2[j] = new Random().nextInt(nFS + 1);
							}
						}
					}
					dFitnLear = getFitnessLearning(offs1, dataTrain, dFuzzySets, nAttr, nLabe, sizeRule);
					offs1[sizeRule - 1] = dFitnLear;
					dFitnLear = getFitnessLearning(offs2, dataTrain, dFuzzySets, nAttr, nLabe, sizeRule);
					offs2[sizeRule - 1] = dFitnLear;
					popuLear_P_O.add(offs1);
					popuLear_P_O.add(offs2);
				}

				//System.out.println("*** New Population ***");
				double miniFitnLear;
				int iMiniFitnLear;
				while (popuLear_P_O.size() != ps) {
					miniFitnLear = Double.MAX_VALUE;
					iMiniFitnLear = -1;
					for(int i = 0; i < popuLear_P_O.size(); i++) {
						dFitnLear = popuLear_P_O.get(i)[sizeRule - 1];
						if(dFitnLear <= miniFitnLear) {
							miniFitnLear = dFitnLear;
							iMiniFitnLear = i;
						}
					}
					popuLear_P_O.remove(iMiniFitnLear);
				}
			}
			
			double maxiFintLear = -1 * Double.MAX_VALUE;
			int iMaxiFitnLear = -1;
			for (int i = 0; i < ps; i++) {
				dFitnLear = popuLear_P_O.get(i)[sizeRule - 1];
				if(dFitnLear >= maxiFintLear) {
					maxiFintLear = dFitnLear;
					iMaxiFitnLear = i;
				}
			}
			
			System.out.println("\tMaximum Fitness: " + popuLear_P_O.get(iMaxiFitnLear)[sizeRule - 1]);
			preliminaryRB.add(new CChrom(popuLear_P_O.get(iMaxiFitnLear)));
			System.out.print("\tFuzzy rule added to preliminary rule base: ");
			showRule(popuLear_P_O.get(iMaxiFitnLear));
			
			int iRuleClas;
			int iPosiClas = sizeRule - 2;
			
			for (int ex = 0; ex < dataTrain.length; ex++) {
				double df_fr = Double.MAX_VALUE;
				for(int i = 0; i < nAttr; i++) {
					dMembDegr = getMembershipGrade(dFuzzySets, dataTrain[ex][i], i, (int) popuLear_P_O.get(iMaxiFitnLear)[i]);
					if (dMembDegr == 0.0) {
						df_fr = 0.0;
						break;
					}
					if (dMembDegr < df_fr) {
						df_fr = dMembDegr;
					}
				}
				iRuleClas = (int) popuLear_P_O.get(iMaxiFitnLear)[iPosiClas];
				if(dataTrain[ex][nAttr + iRuleClas] > 0.0) {
					dataTrain[ex][nAttr + iRuleClas] = dataTrain[ex][nAttr + iRuleClas] - dataTrain[ex][nAttr + iRuleClas] * df_fr;
				}
			}
			popuLear_P_O.clear();
		}
		
		//******* Tuning Preliminary Rule Base ************************
		int sizeChroTunn = preliminaryRB.size() + 2;
		double[][] popuTuni_P_O = new double[pst*2][sizeChroTunn];
		double[][] popuTuni_P_O_aux = new double[pst][sizeChroTunn];
		int iMini1 = -1;
		int iMini2 = -1;
		int iMiniAux;
		double dP1, dP2, dBeta, dRand;
		for(int iC = 0; iC < sizeChroTunn - 1; iC++) {
			popuTuni_P_O[0][iC] = 1.0;
		}
		popuTuni_P_O[0][sizeChroTunn - 2] = 0.0;
		for(int iI = 1; iI < pst; iI++) {
			for(int iC = 0; iC < sizeChroTunn - 1; iC++) {
				popuTuni_P_O[iI][iC] = new Random().nextDouble();
			}
		}
		for(int iI = 0; iI < pst; iI++) {
			dHammTrain = getHammLoss(preliminaryRB, dataTrain, aRealClasTrain, dFuzzySets, nAttr, nLabe, sizeRule, popuTuni_P_O[iI], 1);
			popuTuni_P_O[iI][sizeChroTunn - 1] = dHammTrain;
		}
		
		System.out.println("Genetic tunning of preliminary rule base");
		for(int iT = 0; iT < nit; iT++) {
			System.out.println("\tIteration of genetic algorithm " + (iT + 1) + " of " + nit); 
			for (int iI = pst; iI < pst * 2; iI++) {
				//System.out.println("Iteration " + (iT + 1) + ": Individual " + (iI - iAG_Size));
				for(int t = 0; t < tst; t++) {
					i1 = new Random().nextInt(pst);
					i2 = new Random().nextInt(pst);
					while (i1 == i2) {
						i2 = new Random().nextInt(pst);
					}
					dAccu1 = popuTuni_P_O[i1][sizeChroTunn - 1];
					dAccu2 = popuTuni_P_O[i2][sizeChroTunn - 1];
					if (dAccu1 < dAccu2) {
						iMiniAux = i1;
					} else if (dAccu2 < dAccu1) {
						iMiniAux = i2;
					} else {
						if (new Random().nextDouble() < 0.5) {
							iMiniAux = i1;
						} else {
							iMiniAux = i2;
						}
					}
					if(t == 0) {
						iMini1 = iMiniAux;
					}
					else {
						iMini2 = iMiniAux;
					}
				}
				// Crossover
				for (int j = 0; j < sizeChroTunn - 1; j++) {
					dP1 = popuTuni_P_O[iMini1][j];
					dP2 = popuTuni_P_O[iMini2][j];
					dBeta = GetRandomNumber(dAlfa * -1, 1 + dAlfa);
					popuTuni_P_O[iI][j] =  dP1 + dBeta*(dP2 - dP1);
					if(popuTuni_P_O[iI][j] > 1.0) {
						popuTuni_P_O[iI][j] = 1.0;
					}
					if(popuTuni_P_O[iI][j] < 0.0) {
						popuTuni_P_O[iI][j] = 0.0;
					}
				}
				dRand = new Random().nextDouble();
				if(dRand < mpt) {
					for(int j = 0; j < sizeChroTunn - 2; j++) {
						if(new Random().nextDouble() < mpt_wth) {
							popuTuni_P_O[iI][j] =  new Random().nextDouble();
						}
					}
				}
				dHammTrain = getHammLoss(preliminaryRB, dataTrain, aRealClasTrain, dFuzzySets, nAttr, nLabe, sizeRule, popuTuni_P_O[iI], 1);
				popuTuni_P_O[iI][sizeChroTunn - 1] = dHammTrain;
			}
			//Elistim
			double dAccuMiniTrain;
			int iIndeMiniTrain;
			double dAccuRuleTrain;
			for(int iP = 0; iP < pst; iP++) {
				dAccuMiniTrain = Double.MAX_VALUE;
				iIndeMiniTrain = -1;
				for(int i = 0; i < pst * 2; i++) {
					dAccuRuleTrain = popuTuni_P_O[i][sizeChroTunn - 1];
					if(dAccuRuleTrain < dAccuMiniTrain) {
						dAccuMiniTrain = dAccuRuleTrain;
						iIndeMiniTrain = i;
					}
				}
				for(int j = 0; j < sizeChroTunn; j++) {
					popuTuni_P_O_aux[iP][j] = popuTuni_P_O[iIndeMiniTrain][j];
				}
				popuTuni_P_O[iIndeMiniTrain][sizeChroTunn - 1] = Double.MAX_VALUE;
			}
			
			//New Population
			for(int i = 0; i < pst; i++) {
				for(int j = 0; j < sizeChroTunn; j++) {
					popuTuni_P_O[i][j] = popuTuni_P_O_aux[i][j];
				}
			}
			//Print
			dHammTrain = popuTuni_P_O[0][sizeChroTunn - 1];
			System.out.println("\tMinimum hamming loss: " + dHammTrain);
		}
		
		System.out.println("*** Final Rule Base ******");
		showRB(preliminaryRB, popuTuni_P_O[0], nAttr, nLabe);
		
		System.out.println("*** Final Results ******");
		
		double dClasAccuTrain = getClasAccu(preliminaryRB, dataTrain, aRealClasTrain, dFuzzySets, nAttr, nLabe, sizeRule, popuTuni_P_O[0], 2);
		double dClasAccuTest = getClasAccu(preliminaryRB, dataTest, aRealClasTest, dFuzzySets, nAttr, nLabe, sizeRule, popuTuni_P_O[0], 2);
		System.out.println("Classification Error Train: " + dClasAccuTrain);
		System.out.println("Classification Error Test: " + dClasAccuTest);
		
		dHammTrain = getHammLoss(preliminaryRB, dataTrain, aRealClasTrain, dFuzzySets, nAttr, nLabe, sizeRule, popuTuni_P_O[0], 2);
		dHammTest = getHammLoss(preliminaryRB, dataTest, aRealClasTest, dFuzzySets, nAttr, nLabe, sizeRule, popuTuni_P_O[0], 2);
		System.out.println("Hamming Loss Train: " + dHammTrain);
		System.out.println("Hamming Loss Test: " + dHammTest);
		
		double dAccuTrain = getAccu(preliminaryRB, dataTrain, aRealClasTrain, dFuzzySets, nAttr, nLabe, sizeRule, popuTuni_P_O[0], 2);
		double dAccuTest = getAccu(preliminaryRB, dataTest, aRealClasTest, dFuzzySets, nAttr, nLabe, sizeRule, popuTuni_P_O[0], 2);
		System.out.println("Accuracy Train: " + dAccuTrain);
		System.out.println("Accurayc Test: " + dAccuTest);
		
		double dPrecTrain = getPrec(preliminaryRB, dataTrain, aRealClasTrain, dFuzzySets, nAttr, nLabe, sizeRule, popuTuni_P_O[0], 2);
		double dPrecTest = getPrec(preliminaryRB, dataTest, aRealClasTest, dFuzzySets, nAttr, nLabe, sizeRule, popuTuni_P_O[0], 2);
		System.out.println("Precision Train: " + dPrecTrain);
		System.out.println("Precision Test: " + dPrecTest);
		
		double dRecaTrain = getReca(preliminaryRB, dataTrain, aRealClasTrain, dFuzzySets, nAttr, nLabe, sizeRule, popuTuni_P_O[0], 2);
		double dRecaTest = getReca(preliminaryRB, dataTest, aRealClasTest, dFuzzySets, nAttr, nLabe, sizeRule, popuTuni_P_O[0], 2);
		System.out.println("Recall Train: " + dRecaTrain);
		System.out.println("Recall Test: " + dRecaTest);
		
		double dF1Train = getF1(preliminaryRB, dataTrain, aRealClasTrain, dFuzzySets, nAttr, nLabe, sizeRule, popuTuni_P_O[0], 2);
		double dF1Test = getF1(preliminaryRB, dataTest, aRealClasTest, dFuzzySets, nAttr, nLabe, sizeRule, popuTuni_P_O[0], 2);
		System.out.println("F1-measure Train: " + dF1Train);
		System.out.println("F1-measure Test: " + dF1Test);

		System.out.println("*** END ***");
	}

	private static void showRule(double[] rule) {
		for (int i = 0; i < rule.length - 1; i++) {
			if(i != rule.length - 2) {
				if(rule[i] != 0.0) {
					System.out.print("X" + (i + 1) + " IS " + getLingTerm((int) rule[i]));
					System.out.print(" ");
				}
			}
			else {
				System.out.println("THEN CLASS IS " + rule[i]);
			}
		}
	}

	private static void saveRB(String sDS, int sizeFinalRB, int nFuzzySets, int sizePopu, int iIter,
			ArrayList<CChrom> finalRB)  throws IOException {
		DateTimeFormatter dtf = DateTimeFormatter.ofPattern("yyyy_MM_dd_HH_mm_ss_SSS");  
		LocalDateTime now = LocalDateTime.now();  
		String sFileTrain = "RB_" + sDS + "_" + sizeFinalRB + "_" + nFuzzySets + "_" + sizePopu + "_" + iIter + "_" + dtf.format(now) + ".txt";
		BufferedWriter fRB = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(new File(sFileTrain))));
		double[] dRule;
		for (int i = 0; i < finalRB.size(); i++) {
			dRule = finalRB.get(i).rule;
			for(int c = 0; c < dRule.length - 1; c++) {
				if(c != dRule.length - 2) {
					fRB.write(String.valueOf(dRule[c]) + " ");
				}
				else {
					fRB.write(String.valueOf(dRule[c]));
				}
			}
			fRB.newLine();
		}
		fRB.close();
	}

	private static double getFitnessLearning(double[] rule, double[][] Ddata, double[][] dFuzzySets, int Iinputs, int Ioutputs, int sizeRule) {
		double dMembDegr;
		double dDisparo = 0;
		double dClas = 0.0;
		double dOtherClas = 0.0;
		int iPosiClas = sizeRule - 2;
		int iRuleClas;
		boolean bAllDontCare = true;
		
		for(int i = 0; i < Iinputs; i++) {
			if(rule[i] != 0.0) {
				bAllDontCare = false;
				break;
			}
		}
		if(bAllDontCare) {
			return Double.MAX_VALUE * -1;
		}
		for (int ex = 0; ex < Ddata.length; ex++) {
			dDisparo = Double.MAX_VALUE;
			for(int i = 0; i < Iinputs; i++) {
				dMembDegr = getMembershipGrade(dFuzzySets, Ddata[ex][i], i, (int) rule[i]);
				if (dMembDegr == 0.0) {
					dDisparo = 0.0;
					break;
				}
				if (dMembDegr < dDisparo) {
					dDisparo = dMembDegr;
				}
			}
			iRuleClas = (int) rule[iPosiClas];
			if(Ddata[ex][Iinputs + iRuleClas] > 0.0) {
				dClas = dClas + dDisparo * Ddata[ex][Iinputs + iRuleClas];
			}
			else {
				dOtherClas = dOtherClas + dDisparo;
			}
		}
		return dClas - dOtherClas;
	}

	private static double getErro(ArrayList<CChrom> RB, double[][] Ddata, int[][] aRealClas, double[][] dFuzzySets, int Iinputs, int Ioutputs,
			int sizeRule, int iTrain) {
		int iErro = 0;
		double dMiniMembDegr;
		double dMembDegr;
		int[] aPredClass = new int[Ioutputs];
		double[] dMaxiDisparo = new double[Ioutputs];
		boolean bEqua = true;
		double[] dRule;
		int clasRule;

		for (int r = 0; r < Ddata.length; r++) {
			for (int i = 0; i < Ioutputs; i++) {
				aPredClass[i] = 0;
				dMaxiDisparo[i] = -1 * Double.MAX_VALUE;
			}
			for (int i = 0; i < RB.size(); i++) {
				dRule = RB.get(i).rule;
				dMiniMembDegr = 1.0;
				for(int c_rule = 0; c_rule < sizeRule - 2; c_rule++) {
					dMembDegr = getMembershipGrade(dFuzzySets, Ddata[r][c_rule], c_rule, (int) dRule[c_rule]);
					if (dMembDegr == 0.0) {
						dMiniMembDegr = 0.0;
						break;
					}
					if (dMembDegr < dMiniMembDegr) {
						dMiniMembDegr = dMembDegr;
					}
				}
				clasRule = (int) dRule[sizeRule - 2];
				if (dMiniMembDegr > dMaxiDisparo[clasRule]) {
					dMaxiDisparo[clasRule] = dMiniMembDegr;
				}
			}
			for (int i = 0; i < Ioutputs; i++) {
				if(dMaxiDisparo[i] > 0) {
					aPredClass[i] = 1;
				}
			}
			bEqua = true;
			for (int i = 0; i < Ioutputs; i++) {
				if (aPredClass[i] != aRealClas[r][i]) {
					bEqua = false;
					break;
				}
			}
			if(iTrain == 1) {
				//showResult(aPredClass, aRealClas[r], bEqua);
				//showResultProm(dMaxiDisparo, aPredClass, aRealClas[r], 0.0, bEqua);
			}
			if (!bEqua) {
				iErro++;
			}
		}
		return (double) iErro / Ddata.length;
	}

	private static double getHammLoss(ArrayList<CChrom> RB, double[][] Ddata, int[][] aRealClas, double[][] dFuzzySets, int Iinputs, int Ioutputs,
			int sizeRule, double[] dTH, int iTrain) {
		double dMiniMembDegr;
		double dMembDegr;
		double[] aPredClass = new double[Ioutputs];
		double[] dMaxiDisparo = new double[Ioutputs];
		double[] dCantDisparo = new double[Ioutputs];
		double[] dRule;
		int clasRule;
		double dTotal = 0.0;
		double dT;
		double hl = 0.0;
		double hl_total = 0.0;

		for (int r = 0; r < Ddata.length; r++) {
			dTotal = 0.0;
			for (int i = 0; i < Ioutputs; i++) {
				aPredClass[i] = 0.0;
				dMaxiDisparo[i] = 0.0;
				dCantDisparo[i] = 0.0;
			}
			for (int i = 0; i < RB.size(); i++) {
				dRule = RB.get(i).rule;
				dMiniMembDegr = 1.0;
				dMembDegr = 0.0;
				for(int c_rule = 0; c_rule < sizeRule - 2; c_rule++) {
					dMembDegr = getMembershipGrade(dFuzzySets, Ddata[r][c_rule], c_rule, (int) dRule[c_rule]);
					if (dMembDegr == 0.0) {
						dMiniMembDegr = 0.0;
						break;
					}
					if (dMembDegr < dMiniMembDegr) {
						dMiniMembDegr = dMembDegr;
					}
				}
				dMiniMembDegr = dMiniMembDegr * dTH[i];
				clasRule = (int) dRule[sizeRule - 2];
				if(dMiniMembDegr > 0.0) {
					dMaxiDisparo[clasRule] = dMaxiDisparo[clasRule] + dMiniMembDegr;
				}
			}
			dT = dTH[dTH.length - 2];
			for (int i = 0; i < Ioutputs; i++) {
				dTotal = dTotal + dMaxiDisparo[i];
			}
			for (int i = 0; i < Ioutputs; i++) {
				if(dTotal > 0.0) {
					dMaxiDisparo[i] = dMaxiDisparo[i] / dTotal;
				}
				if(dMaxiDisparo[i] > dT) {
					aPredClass[i] = 1;
				}
			}
			for (int i = 0; i < Ioutputs; i++) {
				hl_total++;
				if (aPredClass[i] != aRealClas[r][i]) {
					hl++;
				}
			}
			if(iTrain == 2) {
				//showResult(aPredClass, aRealClas[r], bEqua);
				//showResultProm(dMaxiDisparo, aPredClass, aRealClas[r], dT, bEqua);
			}
		}
		return hl / hl_total;
	}
	
	private static double getAccu(ArrayList<CChrom> RB, double[][] Ddata, int[][] aRealClas, double[][] dFuzzySets, int Iinputs, int Ioutputs,
			int sizeRule, double[] dTH, int iTrain) {
		double dMiniMembDegr;
		double dMembDegr;
		double[] aPredClass = new double[Ioutputs];
		double[] dMaxiDisparo = new double[Ioutputs];
		double[] dCantDisparo = new double[Ioutputs];
		double[] dRule;
		int clasRule;
		double dTotal = 0.0;
		double dT;
		double dAccuNume = 0.0;
		double dAccuDeno = 0.0;

		for (int r = 0; r < Ddata.length; r++) {
			dTotal = 0.0;
			for (int i = 0; i < Ioutputs; i++) {
				aPredClass[i] = 0.0;
				dMaxiDisparo[i] = 0.0;
				dCantDisparo[i] = 0.0;
			}
			for (int i = 0; i < RB.size(); i++) {
				dRule = RB.get(i).rule;
				dMiniMembDegr = 1.0;
				dMembDegr = 0.0;
				for(int c_rule = 0; c_rule < sizeRule - 2; c_rule++) {
					dMembDegr = getMembershipGrade(dFuzzySets, Ddata[r][c_rule], c_rule, (int) dRule[c_rule]);
					if (dMembDegr == 0.0) {
						dMiniMembDegr = 0.0;
						break;
					}
					if (dMembDegr < dMiniMembDegr) {
						dMiniMembDegr = dMembDegr;
					}
				}
				dMiniMembDegr = dMiniMembDegr * dTH[i];
				clasRule = (int) dRule[sizeRule - 2];
				if(dMiniMembDegr > 0.0) {
					dMaxiDisparo[clasRule] = dMaxiDisparo[clasRule] + dMiniMembDegr; 
				}
			}
			dT = dTH[dTH.length - 2];
			for (int i = 0; i < Ioutputs; i++) {
				dTotal = dTotal + dMaxiDisparo[i];
			}
			for (int i = 0; i < Ioutputs; i++) {
				if(dTotal > 0.0) {
					dMaxiDisparo[i] = dMaxiDisparo[i] / dTotal;
				}
				if(dMaxiDisparo[i] > dT) {
					aPredClass[i] = 1;
				}
			}
			for (int i = 0; i < Ioutputs; i++) {
				if (aRealClas[r][i] == 1.0 &&  aPredClass[i] == 1.0) {
					dAccuNume++;
				}
				if (aRealClas[r][i] == 1.0 || aPredClass[i]  == 1.0) {
					dAccuDeno++;
				}
			}
			if(iTrain == 2) {
				//showResult(aPredClass, aRealClas[r], bEqua);
				//showResultProm(dMaxiDisparo, aPredClass, aRealClas[r], dT, bEqua);
			}
		}
		return dAccuNume / dAccuDeno;
	}
	
	private static double getPrec(ArrayList<CChrom> RB, double[][] Ddata, int[][] aRealClas, double[][] dFuzzySets, int Iinputs, int Ioutputs,
			int sizeRule, double[] dTH, int iTrain) {
		double dMiniMembDegr;
		double dMembDegr;
		double[] aPredClass = new double[Ioutputs];
		double[] dMaxiDisparo = new double[Ioutputs];
		double[] dCantDisparo = new double[Ioutputs];
		double[] dRule;
		int clasRule;
		double dTotal = 0.0;
		double dT;
		double dPrecNume = 0.0;
		double dPrecDeno = 0.0;

		for (int r = 0; r < Ddata.length; r++) {
			dTotal = 0.0;
			for (int i = 0; i < Ioutputs; i++) {
				aPredClass[i] = 0.0;
				dMaxiDisparo[i] = 0.0;
				dCantDisparo[i] = 0.0;
			}
			for (int i = 0; i < RB.size(); i++) {
				dRule = RB.get(i).rule;
				dMiniMembDegr = 1.0;
				dMembDegr = 0.0;
				for(int c_rule = 0; c_rule < sizeRule - 2; c_rule++) {
					dMembDegr = getMembershipGrade(dFuzzySets, Ddata[r][c_rule], c_rule, (int) dRule[c_rule]);
					if (dMembDegr == 0.0) {
						dMiniMembDegr = 0.0;
						break;
					}
					if (dMembDegr < dMiniMembDegr) {
						dMiniMembDegr = dMembDegr;
					}
				}
				dMiniMembDegr = dMiniMembDegr * dTH[i];
				clasRule = (int) dRule[sizeRule - 2];
				if(dMiniMembDegr > 0.0) {
					dMaxiDisparo[clasRule] = dMaxiDisparo[clasRule] + dMiniMembDegr; 
				}
			}
			dT = dTH[dTH.length - 2];
			for (int i = 0; i < Ioutputs; i++) {
				dTotal = dTotal + dMaxiDisparo[i];
			}
			for (int i = 0; i < Ioutputs; i++) {
				if(dTotal > 0.0) {
					dMaxiDisparo[i] = dMaxiDisparo[i] / dTotal;
				}
				if(dMaxiDisparo[i] > dT) {
					aPredClass[i] = 1;
				}
			}
			for (int i = 0; i < Ioutputs; i++) {
				if (aRealClas[r][i] == 1.0 &&  aPredClass[i] == 1.0) {
					dPrecNume++;
				}
				if (aPredClass[i]  == 1.0) {
					dPrecDeno++;
				}
			}
			if(iTrain == 2) {
				//showResult(aPredClass, aRealClas[r], bEqua);
				//showResultProm(dMaxiDisparo, aPredClass, aRealClas[r], dT, bEqua);
			}
		}
		return dPrecNume / dPrecDeno;
	}
	
	private static double getReca(ArrayList<CChrom> RB, double[][] Ddata, int[][] aRealClas, double[][] dFuzzySets, int Iinputs, int Ioutputs,
			int sizeRule, double[] dTH, int iTrain) {
		double dMiniMembDegr;
		double dMembDegr;
		double[] aPredClass = new double[Ioutputs];
		double[] dMaxiDisparo = new double[Ioutputs];
		double[] dCantDisparo = new double[Ioutputs];
		double[] dRule;
		int clasRule;
		double dTotal = 0.0;
		double dT;
		double dRecaNume = 0.0;
		double dRecaDeno = 0.0;

		for (int r = 0; r < Ddata.length; r++) {
			dTotal = 0.0;
			for (int i = 0; i < Ioutputs; i++) {
				aPredClass[i] = 0.0;
				dMaxiDisparo[i] = 0.0;
				dCantDisparo[i] = 0.0;
			}
			for (int i = 0; i < RB.size(); i++) {
				dRule = RB.get(i).rule;
				dMiniMembDegr = 1.0;
				dMembDegr = 0.0;
				for(int c_rule = 0; c_rule < sizeRule - 2; c_rule++) {
					dMembDegr = getMembershipGrade(dFuzzySets, Ddata[r][c_rule], c_rule, (int) dRule[c_rule]);
					if (dMembDegr == 0.0) {
						dMiniMembDegr = 0.0;
						break;
					}
					if (dMembDegr < dMiniMembDegr) {
						dMiniMembDegr = dMembDegr;
					}
				}
				dMiniMembDegr = dMiniMembDegr * dTH[i];
				clasRule = (int) dRule[sizeRule - 2];
				if(dMiniMembDegr > 0.0) {
					dMaxiDisparo[clasRule] = dMaxiDisparo[clasRule] + dMiniMembDegr; 
				}
			}
			dT = dTH[dTH.length - 2];
			for (int i = 0; i < Ioutputs; i++) {
				dTotal = dTotal + dMaxiDisparo[i];
			}
			for (int i = 0; i < Ioutputs; i++) {
				if(dTotal > 0.0) {
					dMaxiDisparo[i] = dMaxiDisparo[i] / dTotal;
				}
				if(dMaxiDisparo[i] > dT) {
					aPredClass[i] = 1;
				}
			}
			for (int i = 0; i < Ioutputs; i++) {
				if (aRealClas[r][i] == 1.0 &&  aPredClass[i] == 1.0) {
					dRecaNume++;
				}
				if (aRealClas[r][i] == 1.0) {
					dRecaDeno++;
				}
			}
			if(iTrain == 2) {
				//showResult(aPredClass, aRealClas[r], bEqua);
				//showResultProm(dMaxiDisparo, aPredClass, aRealClas[r], dT, bEqua);
			}
		}
		return dRecaNume / dRecaDeno;
	}
	
	private static double getF1(ArrayList<CChrom> RB, double[][] Ddata, int[][] aRealClas, double[][] dFuzzySets, int Iinputs, int Ioutputs,
			int sizeRule, double[] dTH, int iTrain) {
		double dMiniMembDegr;
		double dMembDegr;
		double[] aPredClass = new double[Ioutputs];
		double[] dMaxiDisparo = new double[Ioutputs];
		double[] dCantDisparo = new double[Ioutputs];
		double[] dRule;
		int clasRule;
		double dTotal = 0.0;
		double dT;
		double dF1Nume = 0.0;
		double dF1Deno = 0.0;

		for (int r = 0; r < Ddata.length; r++) {
			dTotal = 0.0;
			for (int i = 0; i < Ioutputs; i++) {
				aPredClass[i] = 0.0;
				dMaxiDisparo[i] = 0.0;
				dCantDisparo[i] = 0.0;
			}
			for (int i = 0; i < RB.size(); i++) {
				dRule = RB.get(i).rule;
				dMiniMembDegr = 1.0;
				dMembDegr = 0.0;
				for(int c_rule = 0; c_rule < sizeRule - 2; c_rule++) {
					dMembDegr = getMembershipGrade(dFuzzySets, Ddata[r][c_rule], c_rule, (int) dRule[c_rule]);
					if (dMembDegr == 0.0) {
						dMiniMembDegr = 0.0;
						break;
					}
					if (dMembDegr < dMiniMembDegr) {
						dMiniMembDegr = dMembDegr;
					}
				}
				dMiniMembDegr = dMiniMembDegr * dTH[i];
				clasRule = (int) dRule[sizeRule - 2];
				if(dMiniMembDegr > 0.0) {
					dMaxiDisparo[clasRule] = dMaxiDisparo[clasRule] + dMiniMembDegr; 
				}
			}
			dT = dTH[dTH.length - 2];
			for (int i = 0; i < Ioutputs; i++) {
				dTotal = dTotal + dMaxiDisparo[i];
			}
			for (int i = 0; i < Ioutputs; i++) {
				if(dTotal > 0.0) {
					dMaxiDisparo[i] = dMaxiDisparo[i] / dTotal;
				}
				if(dMaxiDisparo[i] > dT) {
					aPredClass[i] = 1;
				}
			}
			for (int i = 0; i < Ioutputs; i++) {
				if (aRealClas[r][i] == 1.0 &&  aPredClass[i] == 1.0) {
					dF1Nume++;
				}
				if (aRealClas[r][i] == 1.0) {
					dF1Deno++;
				}
				if (aPredClass[i] == 1.0) {
					dF1Deno++;
				}
			}
			if(iTrain == 2) {
				//showResult(aPredClass, aRealClas[r], bEqua);
				//showResultProm(dMaxiDisparo, aPredClass, aRealClas[r], dT, bEqua);
			}
		}
		return 2 * dF1Nume / dF1Deno;
	}

	private static double getClasAccu(ArrayList<CChrom> RB, double[][] Ddata, int[][] aRealClas, double[][] dFuzzySets, int Iinputs, int Ioutputs,
			int sizeRule, double[] dTH, int iTrain) {
		int iErro = 0;
		double dMiniMembDegr;
		double dMembDegr;
		int[] aPredClass = new int[Ioutputs];
		double[] dMaxiDisparo = new double[Ioutputs];
		double[] dCantDisparo = new double[Ioutputs];
		boolean bEqua = true;
		double[] dRule;
		int clasRule;
		double dTotal = 0.0;
		int iTotal = 0;
		double dT;

		for (int r = 0; r < Ddata.length; r++) {
			dTotal = 0.0;
			iTotal = 0;
			for (int i = 0; i < Ioutputs; i++) {
				aPredClass[i] = 0;
				dMaxiDisparo[i] = 0.0;
				dCantDisparo[i] = 0.0;
			}
			for (int i = 0; i < RB.size(); i++) {
				dRule = RB.get(i).rule;
				dMiniMembDegr = 1.0;
				dMembDegr = 0.0;
				for(int c_rule = 0; c_rule < sizeRule - 2; c_rule++) {
					dMembDegr = getMembershipGrade(dFuzzySets, Ddata[r][c_rule], c_rule, (int) dRule[c_rule]);
					if (dMembDegr == 0.0) {
						dMiniMembDegr = 0.0;
						break;
					}
					if (dMembDegr < dMiniMembDegr) {
						dMiniMembDegr = dMembDegr;
					}
				}
				dMiniMembDegr = dMiniMembDegr * dTH[i];
				clasRule = (int) dRule[sizeRule - 2];
				//if (dMiniMembDegr > dMaxiDisparo[clasRule]) {
				//	dMaxiDisparo[clasRule] = dMiniMembDegr;
				//}
				if(dMiniMembDegr > 0.0) {
					dMaxiDisparo[clasRule] = dMaxiDisparo[clasRule] + dMiniMembDegr;
					//dCantDisparo[clasRule] = dCantDisparo[clasRule] + 1.0; 
				}
			}
			dT = dTH[dTH.length - 2];
			for (int i = 0; i < Ioutputs; i++) {
				dTotal = dTotal + dMaxiDisparo[i];
			}
			for (int i = 0; i < Ioutputs; i++) {
				if(dTotal > 0.0) {
					dMaxiDisparo[i] = dMaxiDisparo[i] / dTotal;
				}
				if(dMaxiDisparo[i] > dT) {
					aPredClass[i] = 1;
				}
			}
			bEqua = true;
			for (int i = 0; i < Ioutputs; i++) {
				if (aPredClass[i] != aRealClas[r][i]) {
					bEqua = false;
					break;
				}
			}
			if(iTrain == 2) {
				//showResult(aPredClass, aRealClas[r], bEqua);
				//showResultProm(dMaxiDisparo, aPredClass, aRealClas[r], dT, bEqua);
			}
			if (!bEqua) {
				iErro++;
			}
		}
		return (double) iErro / Ddata.length;
	}
	
	private static void showResultProm(double[] promClas,  int[] aPredClass, int[] aRealClas, double dTH, boolean bStat) {
		for(int i = 0; i < promClas.length; i++) {
			System.out.print(promClas[i] + " ");
		}
		System.out.print(" - ");
		for(int i = 0; i < aPredClass.length; i++) {
			System.out.print(aPredClass[i]);
		}
		System.out.print(" - ");
		System.out.print(dTH);
		System.out.print(" - ");
		for(int i = 0; i < aRealClas.length; i++) {
			System.out.print(aRealClas[i]);
		}
		if(bStat) {
			System.out.print(" - C");
		}
		else {
			System.out.print(" - E");
		}
		System.out.println();
	}

	private static void showResult(int[] aPredClass, int[] aRealClas, boolean bEqua) {
		for(int i = 0; i < aPredClass.length; i++) {
			System.out.print(aPredClass[i]);
		}
		System.out.print(" - ");
		for(int i = 0; i < aRealClas.length; i++) {
			System.out.print(aRealClas[i]);
		}
		if(bEqua) {
			System.out.println(" - C");
		}
		else {
			System.out.println(" - E");
		}
	}

	private static int getActiCondi(int[] RB, int Iinputs, int Ioutputs, int sizeRule, int sizeRuleBase) {
		int iPosiActiRule = Iinputs + Ioutputs;
		int iActiCondi = 0;
		for (int j = 0; j < sizeRuleBase; j++) {
			for (int c_rule = 0; c_rule < Iinputs; c_rule++) {
				if (RB[(j * sizeRule) + c_rule] != 0 && RB[(j * sizeRule) + iPosiActiRule] == 1) {
					iActiCondi++;
				}
			}
		}
		return iActiCondi;
	}

	private static int getRuleActiveRules(int[] RB, int Iinputs, int Ioutputs, int sizeRule, int sizeRuleBase) {
		int iPosiActiRule = Iinputs + Ioutputs;
		int iActiRule = 0;
		for (int j = 0; j < sizeRuleBase; j++) {
			if (RB[(j * sizeRule) + iPosiActiRule] == 1) {
				iActiRule++;
			}
		}
		return iActiRule;

	}

	private static void showRB(ArrayList<CChrom> finalRB, double[] aWeight, int nAttr, int nLabe) {
		double[] rule;
		for(int iRB = 0; iRB < finalRB.size(); iRB++) {
			rule = finalRB.get(iRB).rule;
			for (int i = 0; i < rule.length - 1; i++) {
				if(i != rule.length - 2) {
					if(rule[i] != 0.0) {
						System.out.print("X" + (i + 1) + " IS " + getLingTerm((int) rule[i]));
						System.out.print(" ");
					}
				}
				else {
					System.out.print("THEN CLASS IS " + rule[i]);
				}
			}
			System.out.println(" (Weight = " + aWeight[iRB] + ")");
		}
	}
	
	private static String getActive(int i) {
		switch (i) {
		case 0:
			return "No Active";
		case 1:
			return "Active";
		default:
			return "Undefined";
		}
	}

	private static String getLingTerm(int lt) {
		String sLingTerm = "?";

		switch (lt) {
		case 0:
			sLingTerm = "Don't Care";
			break;
		case 1:
			sLingTerm = "Very Low";
			break;
		case 2:
			sLingTerm = "Low";
			break;
		case 3:
			sLingTerm = "Medium";
			break;
		case 4:
			sLingTerm = "High";
			break;
		case 5:
			sLingTerm = "Very High";
			break;
		}

		return sLingTerm;
	}

	private static double getMembershipGrade(double[][] dFuzzySets, double x, int iLingTerm, int iFuzzySet) {
		double a;
		double b;
		double c;

		if (iFuzzySet == 0) {
			return 1.0;
		}

		iFuzzySet--;
		a = dFuzzySets[iLingTerm][iFuzzySet * 3];
		b = dFuzzySets[iLingTerm][iFuzzySet * 3 + 1];
		c = dFuzzySets[iLingTerm][iFuzzySet * 3 + 2];

		//System.out.println(a + "\t" + b + "\t" + c);

		if (x <= a) {
			return 0.0;
		}
		if (x > a && x <= b) {
			return (x - a) / (b - a);
		}
		if (x > b && x < c) {
			return (c - x) / (c - b);
		}
		if (x >= c) {
			return 0.0;
		}

		return 0.0;
	}

	private static void showMatrizDouble(double[][] matrix) {
		for (int r = 0; r < matrix.length; r++) {
			for (int c = 0; c < matrix[r].length; c++) {
				System.out.print(matrix[r][c] + "\t");
			}
			System.out.println();
		}
	}

	private static void showMatrizInteger(int[][] matrix) {
		for (int r = 0; r < matrix.length; r++) {
			for (int c = 0; c < matrix[r].length; c++) {
				System.out.print(matrix[r][c] + "\t");
			}
			System.out.println();
		}
	}
	
	public static double GetRandomNumber(double minimum, double maximum)
	{ 
	    Random random = new Random();
	    return random.nextDouble() * (maximum - minimum) + minimum;
	}

}