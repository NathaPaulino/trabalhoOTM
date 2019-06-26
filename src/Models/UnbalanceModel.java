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

 *			Y_{j} =  \begin{Bmatrix}
				0 & se & (\sum_{a=1}^{n} L_{a} *y_{ia} * g_{a}) + K > 0 \\ 
				1 & se & (\sum_{a=1}^{n} L_{a} *y_{ia} * g_{a}) + K \leq 0
			\end{Bmatrix} \forall j
			
 * 		Restrições	
 			g_{n} \in \{0,1\} \: \: \forall n \: n \in 1...n
  			X_{i} \in \{0,1\} \: \: \forall i \: i \in 1...PCR
 			Y_{j} \in \{0,1\} \: \: \forall j \: j \in 1...RD
 * 
 * 		Modelo de Otimização - II: 
 * 		
 * 		\[\max Z(p_{i}^{+},p_{i}^{-},p_{j}^{+},p_{j}^{-},g_{a}) = (\sum_{i=1}^{pcr} -(p_{i}^{+}-p_{i}^{-}))+(\sum_{j=1}^{rd} p_{j}^{+} - p_{j}^{-}))\]
		
		Sujeito a:
			
			\[p_{j}^{+} - p_{j}^{-} = (\sum_{a=1}^{n}g_{a}*L_{a}*y_{ja})+K, \, \, \, \forall j, \,\, j \in \[1,rd\]\]
			\[p_{i}^{+} - p_{i}^{-} = (\sum_{a=1}^{n}g_{a}*L_{a}*x_{ia})+K, \, \, \, \forall i, \,\, i \in \[1,pcr\]\]
			\[p_{i}^{+},p_{i}^{-},p_{j}^{+},p_{j}^{-} \geqslant 0, \forall i,\forall j, i \in \[1,pcr\] \, , j \in [1,rd]\]
			g_{a} \in \{0,1\}, \, \, \forall a , \, \, a \in \[1,n\]
 */

public class UnbalanceModel {
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
				System.out.println("Solução Encontrada! - Amostra Desbalanceada!");
				System.out.println("Objetivo = "+ cplex.getObjValue());
				for(a=0;a<n;a++) {
					System.out.println("G["+(a+1)+"]: "+cplex.getValue(g[a]));
				}
				for(a=0;a<constraints.size();a++) {
					System.out.println("Slack constraints: "+cplex.getSlack(constraints.get(a)));
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
	     * "/home/np/eclipse-workspace/TrabalhoOTM/csvFiles/UnbalancedData/UnbalanceCoeff.csv" --> L
	     * "/home/np/eclipse-workspace/TrabalhoOTM/csvFiles/UnbalancedData/UnbalanceCoeffK.csv" --> K
	     * "/home/np/eclipse-workspace/TrabalhoOTM/csvFiles/UnbalancedData/UnbalanceTestPCR.csv" --> TestPCR
	     * "/home/np/eclipse-workspace/TrabalhoOTM/csvFiles/UnbalancedData/UnbalanceTestRD.csv" --> TestRD
	     */ 
	/*	File csvFile  = new File("/home/np/eclipse-workspace/TrabalhoOTM/csvFiles/UnbalancedData/UnbalanceCoeff.csv");
		File csvFile2 = new File("/home/np/eclipse-workspace/TrabalhoOTM/csvFiles/UnbalancedData/UnbalanceCoeffK.csv");
		File csvFile3 = new File("/home/np/eclipse-workspace/TrabalhoOTM/csvFiles/UnbalancedData/UnbalanceTestPCR.csv");
		File csvFile4 = new File("/home/np/eclipse-workspace/TrabalhoOTM/csvFiles/UnbalancedData/UnbalanceTestRD.csv");
	*/	File csvFile  = new File("/home/np/eclipse-workspace/TrabalhoOTM/csvFiles/UnbalancedData/UnbalanceCoeff.csv");
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
			try(BufferedReader csvReader4 = new BufferedReader(new FileReader(csvFile4))){
			   // Parse da matriz de coeficientes relativo aos pacientes RD
			   int i = 0;
			   int j = 0;
			   while ((line = csvReader4.readLine()) != null) {
					String[] values = line.split(",");
					for(j=0;j< values.length;j++) {
						y[i][j] = new Double(values[j]);
					}
					i++;
				}
			}catch (IOException exec) {
				exec.printStackTrace();
			}
		/*String line = "";
	    
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
			}*/
			System.out.println("Valor de k: "+k);
			//Função model(geneNumber,pcr,vectorpcr,rd,vectorRD,paramL,paramK)
			UnbalanceModel.model(n,pcr,x,rd,y,l,k);
		}		
	}
	public static void unbalancedParse() {
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
			try(BufferedReader csvReader4 = new BufferedReader(new FileReader(csvFile4))){
			   // Parse da matriz de coeficientes relativo aos pacientes RD
			   int i = 0;
			   int j = 0;
			   while ((line = csvReader4.readLine()) != null) {
					String[] values = line.split(",");
					for(j=0;j< values.length;j++) {
						y[i][j] = new Double(values[j]);
					}
					i++;
				}
			}catch (IOException exec) {
				exec.printStackTrace();
			}
			UnbalanceModel.model(n, pcr, x, rd, y, l, k);
		}		
	}
}
