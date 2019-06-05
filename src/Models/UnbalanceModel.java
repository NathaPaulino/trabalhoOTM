package Models;

import java.io.*;
import java.util.*;
import ilog.cplex.*;
import ilog.concert.*;

/*
 * 	O DLDA - Diagonal Linear Discriminant Analysis é calculado através da seguinte fórmula:
 *  	DLDA(x) = L*x + K, onde:
 *  		L(i) é peso que o gene i tem na classificação
 *  		x(i,j) é o valor da expressão gênica do gene i no paciente j
 *  		K é uma constante criada pelo próprio modelo.
 *  
 *  Contextualizando para o modelo, a nova equação que delimita o acerto é vinculada por:
 *  	DLDA(x,g) = g*(L*x + K), onde
 *  		L(i)   é peso que o gene i tem na classificação
 *  		x(i,j) é o valor da expressão gênica do gene i no paciente j
 *  		g(i)   é a variável de decisão que identifica se aquele gene é importante na classificação
 *  		K 	   é a constante criada pelo modelo.   
 *	
 *	Modelo de Otimização:
 *		
 *		\max Z(g_{i}) = \sum_{i=1}^{pcr}X_{i} + \sum_{j=1}^{rd}Y_{j}    
 *		Onde:
 *		    X_{i} = \begin{Bmatrix}
				0 & se & (\sum_{a=1}^{n} L_{a} *x_{ia} * g_{a}) + K > 0 \\ 
				1 & se & (\sum_{a=1}^{n} L_{a} *x_{ia} * g_{a}) + K \leq 0
		    	\end{Bmatrix} \forall i

 *			Y_{j} =  	
 *			\begin{Bmatrix}
				0 & se & (\sum_{a=1}^{n} L_{a} *y_{ia} * g_{a}) + K > 0 \\ 
				1 & se & (\sum_{a=1}^{n} L_{a} *y_{ia} * g_{a}) + K \leq 0
			\end{Bmatrix} \forall j
			
 * 		Restrições	
 * 			\sum_{a=1}^{n} g_{a} \leq 81
 * 			\sum_{i=1}^{pcr} X_{i} \leq PCR
 * 			\sum_{j=1}^{rd} Y_{j} \leq RD
 * 			g_{n} \in \{0,1\} \: \: \forall n \: n \in 1...81
 * 			X_{i} \in \{0,1\} \: \: \forall i \: i \in 1...PCR
 *			Y_{j} \in \{0,1\} \: \: \forall j \: j \in 1...RD
 * 
 */

public class UnbalanceModel {
	public static void model(int n,int pcr,double[][] x,int rd,double[][]y,double []L,double K) {
		try{
			
			//Declaração de variáveis,expressões e companhia
			int a,i,j;

			IloCplex cplex = new IloCplex();
			
			IloNumVar []g = new IloNumVar[n];
			IloNumVar []X = new IloNumVar[pcr];
			IloNumVar []Y = new IloNumVar[rd];
			
			IloLinearNumExpr genes		 	 = cplex.linearNumExpr();
			IloLinearNumExpr pcrPacients 	 = cplex.linearNumExpr();
			IloLinearNumExpr rdPacients  	 = cplex.linearNumExpr();
			IloLinearNumExpr pcrEvaluation   = cplex.linearNumExpr();
			IloLinearNumExpr rdEvaluation	 = cplex.linearNumExpr();
			
			for(a=1;a<=n;a++) {
				g[a] = cplex.boolVar();
				genes.addTerm(1.0,g[a]);
			}
			for(i=1;i<=pcr;i++) {
				X[i] = cplex.boolVar();
				pcrPacients.addTerm(1.0, X[i]);
			}
			for(j=1;j<=rd;j++) {
				Y[i] = cplex.boolVar();
				rdPacients.addTerm(1.0, Y[i]);
			} 
			
			for(i=1;i<=pcr;i++) {
				for(a=1;a<=n;a++) {
					pcrEvaluation.addTerm(g[a],L[a]*x[i][a]);	
				}
				pcrEvaluation.addTerm(K,null);
				//
				/*if(pcrEvaluation. 0) {
					X[i]
				}*/
			}
			
			for(j=1;j<=rd;j++) {
				for(a=1;a<=n;a++) {
				   rdEvaluation.addTerm(g[a],L[a]*y[i][a]);
				}
				rdEvaluation.addTerm(K, null);
				/*if(rdEvaluation > 0) {
					Y[i] = 0;
				}*/
			}
			
			//Função Objetivo
			
			IloLinearNumExpr Z = cplex.linearNumExpr();
			
			for(i=1;i<=pcr;i++) {
				Z.addTerm(1.0,X[i]);
			}
			Z.addTerm(rd, null);
			for(j=1;j<=rd;j++) {
				Z.addTerm(-1.0, Y[i]);
			}
			
			//Definindo o tipo da função
			cplex.addMaximize(Z);
			
			// Definindo constraints 
			
			List<IloRange> constraints   = new ArrayList<IloRange>();
			
			constraints.add(cplex.addLe(genes, 81 ));
			constraints.add(cplex.addLe(pcrPacients, pcr));
			constraints.add(cplex.addLe(rdPacients, rd ));
						
			//Se possível resolver o sistema
			if(cplex.solve()) {
				System.out.println("Objetivo = "+ cplex.getObjValue());
				System.out.println("Verificando contas: ");
				for (i=1;i<=pcr;i++) {
					System.out.println("X["+i+"] = "+ cplex.getValue(X[i]));
				}
				for (j=1;j<=rd;j++) {
					System.out.println("Y["+j+"] = "+ cplex.getValue(Y[j]));
				}
				for (a=1;a<=n;a++) {
					System.out.println("g["+a+"] = "+ cplex.getValue(g[a]));
				}
				for (int k=0;k<constraints.size();k++) {
					System.out.println("Dual  constraint: "+(k+1)+" = "+cplex.getDual(constraints.get(i)));
					System.out.println("Slack constraint: "+(k+1)+" = "+cplex.getSlack(constraints.get(i)));
				}
			} else {
				System.out.println("Modelo não resolvido");
			}
			
		} catch (IloException exec) {
			exec.printStackTrace();
		}
		
	}
	
	public static void parseCSV() {
	    /*
	     * "/home/np/eclipse-workspace/TrabalhoOTM/csvFiles/UnbalancedData/UnbalanceCoeff.csv" --> L
	     * "/home/np/eclipse-workspace/TrabalhoOTM/csvFiles/UnbalancedData/UnbalanceCoeffK.csv" --> K
	     * "/home/np/eclipse-workspace/TrabalhoOTM/csvFiles/UnbalancedData/UnbalanceTestPCR.csv" --> TestPCR
	     * "/home/np/eclipse-workspace/TrabalhoOTM/csvFiles/UnbalancedData/UnbalanceTestRD.csv" --> TestRD
	     */ 
		File csvFile  = new File("/home/np/eclipse-workspace/TrabalhoOTM/csvFiles/UnbalancedData/UnbalanceCoeff.csv");
		File csvFile2 = new File("/home/np/eclipse-workspace/TrabalhoOTM/csvFiles/UnbalancedData/UnbalanceCoeffK.csv");
		File csvFile3 = new File("/home/np/eclipse-workspace/TrabalhoOTM/csvFiles/UnbalancedData/UnbalanceTestPCR.csv");
		File csvFile4 = new File("/home/np/eclipse-workspace/TrabalhoOTM/csvFiles/UnbalancedData/UnbalanceTestRD.csv");
	    String line = "";
	    
	    int n = 81;
	    int pcr = 13;
	    int rd = 38;
	    double k=3;
	    double[][] x = new double[pcr][n];
	    double[][] y = new double[rd][n];
	    double [] l = new double [n];
		if(csvFile.isFile() && csvFile2.isFile() && csvFile3.isFile() && csvFile4.isFile()) {
			try(BufferedReader csvReader  = new BufferedReader(new FileReader(csvFile))){ 
				// Parse do coeficiente L do DLDA
				while((line = csvReader.readLine()) != null) {
					String[] values = line.split(",");
					for(int i=0;i <values.length;i++) {
						l[i] = new Double(values[i]);
					}
				}
				System.out.println("Coeficiente L pronto");
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
				System.out.println("Coeficiente K pronto");
			}catch (IOException exec) {
				exec.printStackTrace();
			}
			try(BufferedReader csvReader3 = new BufferedReader(new FileReader(csvFile3))){
			   // Parse da matriz de coeficientes relativo aos pacientes PCR
				int i = 0;
				int j = 0;
				while ((line = csvReader3.readLine()) != null) {
					String[] values = line.split(",");
					System.out.println("Linha: " + i);
					System.out.println("Values length: " + values.length);
					for(j=0;j< values.length;j++) {
						x[i][j] = new Double(values[j]);
					}
					i++;
				}
				System.out.println("Matriz pacientes PCR pronta");
			}catch (IOException exec) {
				exec.printStackTrace();
			}
			try(BufferedReader csvReader4 = new BufferedReader(new FileReader(csvFile4))){
			   // Parse da matriz de coeficientes relativo aos pacientes RD
			   int i = 0;
			   int j = 0;
			   while ((line = csvReader4.readLine()) != null) {
					String[] values = line.split(",");
					System.out.println("Linha: " + i);
					System.out.println("Values length: " + values.length);
					for(j=0;j< values.length;j++) {
						y[i][j] = new Double(values[j]);
					}
					i++;
				}
				System.out.println("Matriz pacientes RD pronta");
			}catch (IOException exec) {
				exec.printStackTrace();
			}
			System.out.println("Valor de k: "+k);
			//Função model(geneNumber,pcr,vectorpcr,rd,vectorRD,paramL,paramK)
			UnbalanceModel.model(n,pcr,x,rd,y,l,k);
		}		
	}
}
