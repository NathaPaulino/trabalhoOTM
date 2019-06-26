package Models;

import java.util.*;
import ilog.cplex.*;
import ilog.concert.*;

public class Example {
	public static void model() {
		try {
			IloCplex cplex = new IloCplex();
			/*
			 *  Criação das variáveis de decisão
			 */
			IloNumVar x = cplex.numVar(0, Double.MAX_VALUE, "x");
			IloNumVar y = cplex.numVar(0, Double.MAX_VALUE, "y");
			
			/*
			 *  Função objetivo
			 */
			 	
			IloLinearNumExpr objetivo = cplex.linearNumExpr();
			objetivo.addTerm(0.12,x);
			objetivo.addTerm(0.15,y);
			//0,12x + 0,15y
			cplex.addMinimize(objetivo);
			
			/*
			 *  Constraints
			 */
			
			List<IloRange> constraints = new ArrayList<IloRange>();
			constraints.add(cplex.addGe(cplex.sum(cplex.prod(60, x),cplex.prod(60, y)),300));
			constraints.add(cplex.addGe(cplex.sum(cplex.prod(12, x),cplex.prod(6, y)),36));
			constraints.add(cplex.addGe(cplex.sum(cplex.prod(10, x),cplex.prod(30, y)),90));
			IloLinearNumExpr num_expr = cplex.linearNumExpr();
			num_expr.addTerm(2,x);
			num_expr.addTerm(-1,y);
			constraints.add(cplex.addEq(num_expr, 0));
			num_expr = cplex.linearNumExpr();
			num_expr.addTerm(1,y);
			num_expr.addTerm(-1,x);
			constraints.add(cplex.addLe(num_expr, 8));
			
			/*
			 * Resolução do modelo e mostrando 
			 * o problema dual com suas constraints.
			 */
			if(cplex.solve()) {
				System.out.println("Objetivo = "+cplex.getObjValue());
				System.out.println("x        = "+cplex.getValue(x));
				System.out.println("y        = "+cplex.getValue(y));
				for (int i=0;i<constraints.size();i++) {
					System.out.println("Dual  constraint: "+(i+1)+" = "+cplex.getDual(constraints.get(i)));
					System.out.println("Slack constraint: "+(i+1)+" = "+cplex.getSlack(constraints.get(i)));
				}
			} else {
				System.out.println("Modelo não resolvido");
			}
			
			cplex.end();
		}catch (IloException exec) {
			exec.printStackTrace();
		}
	}
}
