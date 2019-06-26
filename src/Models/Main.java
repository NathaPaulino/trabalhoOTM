package Models;
import java.util.*;
import java.io.*;

public class Main {

	public static void main(String[] args) {
	/*
	 * O programa funciona invocando a função parser dos arquivos de cada classe.
	  * elisangelamartins@cefetmg.br - Administração
	*/
	/*		  0 	se pi <=0
	 * 	Pi+ = pi	se pi >0
	 * 
	 *		   0  	se pi >=0
	 *  Pi- = -pi	se pi <0
	 *  
	 *		   0    se  pj >0  
	 *  Pj+	=  pj   se  pj <=0
	 *  
	 *         0    se  pj <0
	 *  Pj- =  -pj  se  pj >=0
	 */
	//	Teste.unbalancedParse();
		UnbalanceModel.parseCSV();
		UnderbalanceModel.parseCSV();
		OverbalanceModel.parseCSV();
	}
}
