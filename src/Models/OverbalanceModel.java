package Models;
import java.io.*;
import java.util.*;

import ilog.cplex.*;
import ilog.concert.*;

public class OverbalanceModel {
	public static void model(int n,int pcr,double[][] x,int rd,double[][]y,double []L,double K) {
		try{
			int i,j,a;
			/*
			 *  Instância menor do problema
			 */
			IloCplex cplex = new IloCplex();
			// Variáveis
			IloNumVar []g   = cplex.boolVarArray(n);
			IloNumVar []Pip = cplex.numVarArray(pcr, 0, Double.MAX_VALUE);
			IloNumVar []Pin = cplex.numVarArray(pcr, 0, Double.MAX_VALUE);
			IloNumVar []Pjp = cplex.numVarArray(rd , 0, Double.MAX_VALUE);
			IloNumVar []Pjn = cplex.numVarArray(rd , 0, Double.MAX_VALUE);
			
			// Objetivo
			IloLinearNumExpr objetivo = cplex.linearNumExpr();
			for(i=0;i<pcr;i++) {
				objetivo.addTerm(-1.0, Pip[i]);
				objetivo.addTerm( 1.0, Pin[i]);
			}
			for(j=0;j<rd;j++) {
				objetivo.addTerm( 1.0, Pjp[j]);
				objetivo.addTerm(-1.0, Pjn[j]);
			}
			cplex.addMaximize(objetivo);
			
			// Constraints
			List<IloRange> constraints = new ArrayList<IloRange>();
			
			for(i=0;i<pcr;i++) {
				IloLinearNumExpr difPi = cplex.linearNumExpr();
				difPi.addTerm( 1.0, Pip[i]);
				difPi.addTerm(-1.0, Pin[i]);
				IloLinearNumExpr sumPi = cplex.linearNumExpr();
				for(a=0;a<n;a++) {
					sumPi.addTerm(L[a]*x[i][a], g[a]);
				}
				sumPi.setConstant(K);
				constraints.add((IloRange) cplex.addEq(difPi,sumPi));
			}
			
			for(j=0;j<rd;j++) {
				IloLinearNumExpr difPj = cplex.linearNumExpr();
				difPj.addTerm( 1.0, Pjp[j]);
				difPj.addTerm(-1.0, Pjn[j]);
				IloLinearNumExpr sumPj = cplex.linearNumExpr();
				for(a=0;a<n;a++) {
					sumPj.addTerm(L[a]*y[j][a], g[a]);
				}
				sumPj.setConstant(K);
				constraints.add((IloRange) cplex.addEq(difPj, sumPj));
			}
			if(cplex.solve()) {
				System.out.println("Solução Encontrada! - Amostra Overbalanced");
				System.out.println("Objetivo = "+ cplex.getObjValue());
				for(a=0;a<n;a++) {
					System.out.println("G["+(a+1)+"]: "+cplex.getValue(g[a]));
				}
				for (a=0;a<constraints.size();a++) {
				//	System.out.println("Dual  constraint: "+(c+1)+" = "+cplex.getDual(constraints.get(c)));
					System.out.println("Slack constraint: "+cplex.getSlack(constraints.get(a)));
				}
			} else {
				System.out.println("Modelo não resolvido");
			}
			cplex.end();		
		} catch (IloException exec) {
			exec.printStackTrace();
		}
		
		
	}
	
	public static void parseCSV() {
	    /*
	     * "/home/np/eclipse-workspace/TrabalhoOTM/csvFiles/OverbalancedData/OverbalanceCoeff.csv" --> L
	     * "/home/np/eclipse-workspace/TrabalhoOTM/csvFiles/OverbalancedData/OverbalanceCoeffK.csv" --> K
	     * "/home/np/eclipse-workspace/TrabalhoOTM/csvFiles/OverbalancedData/OverbalanceTestPCR.csv" --> TestPCR
	     * "/home/np/eclipse-workspace/TrabalhoOTM/csvFiles/OverbalancedData/OverbalanceTestRD.csv" --> TestRD
	     */ 
		File csvFile  = new File("/home/np/eclipse-workspace/TrabalhoOTM/csvFiles/OverbalancedData/OverbalanceCoeff.csv");
		File csvFile2 = new File("/home/np/eclipse-workspace/TrabalhoOTM/csvFiles/OverbalancedData/OverbalanceCoeffK.csv");
		File csvFile3 = new File("/home/np/eclipse-workspace/TrabalhoOTM/csvFiles/OverbalancedData/OverbalanceTestPCR.csv");
		File csvFile4 = new File("/home/np/eclipse-workspace/TrabalhoOTM/csvFiles/OverbalancedData/OverbalanceTestRD.csv");
	    String line = "";
	    
	    int n = 81;
	    int pcr = 38;
	    int rd = 38;
	    double k=1.0;
	    
	    double[][] x = new double[pcr][n];
	    double[][] y = new double[rd][n];
	    double [] l = new double [n];
		if(csvFile.isFile() && csvFile2.isFile() && csvFile3.isFile() && csvFile4.isFile()) {
			try(BufferedReader csvReader  = new BufferedReader(new FileReader(csvFile))){ 
				// Parse do coeficiente L do DLDA
				int a=0;
				while((line = csvReader.readLine()) != null) {
					String[] values = line.split(",");
					for(int i=0;i <values.length;i++) {
						l[a] = new Double(values[i]);
						a++;
					}
				}
			}catch (IOException exec) {
				exec.printStackTrace();
			}
			try(BufferedReader csvReader2 = new BufferedReader(new FileReader(csvFile2))){
				// Parse do coeficiente K do DLDA
				while((line = csvReader2.readLine())!= null) {
					String[] values = line.split(",");
					for(int i=0;i<values.length;i++) {
						k = new Double(values[i]);
					}
				}
			}catch (IOException exec) {
				exec.printStackTrace();
			}
			try(BufferedReader csvReader3 = new BufferedReader(new FileReader(csvFile3))){
			   // Parse da matriz de coeficientes relativo aos pacientes PCR
				int i = 0;
				int j = 0;
				while ((line = csvReader3.readLine()) != null) {
					String[] values = line.split(",");
					for(j=0;j< values.length;j++) {
						x[i][j] = new Double(values[j]);
					}
					i++;
				}
			}catch (IOException exec) {
				exec.printStackTrace();
			}
			System.out.println("Valor de k: "+k);
			//Função model(geneNumber,pcr,vectorpcr,rd,vectorRD,paramL,paramK)
			OverbalanceModel.model(n,pcr,x,rd,y,l,k);
		}		
	}

}
