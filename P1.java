
import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Scanner;

class P1
{
	static Scanner input = new Scanner(System.in);
	private static final double TWO_PI = 2.0 * Math.PI;
	private static final double FOUR_PI = 4.0 * Math.PI;
	ArrayList<Double> al = new ArrayList<Double>();
	static String[] splitStr;
	static BufferedReader br =  new BufferedReader (new InputStreamReader(System.in));
	
	public static void main(String args[])
	{	double A [][]= new double[10][10];
		int r= 0, c=0;
		String str1; 
		try {
			while (true) {
				str1 = br.readLine();
				if(str1 == null || str1.isEmpty() || str1.equals("\n") || str1.equals("\r\n")) {
					break;
				}
				str1=str1.trim();
				if (str1.charAt(0) == '#') {
					continue;
				}
				str1=str1.replaceAll("[^-0-9.\\s]", "");
				str1 = str1.replaceAll("\\s+", " ")
						.trim();
				
				splitStr = str1.split(" ");
				int tempCol = 0;
				for (String s : splitStr) {
					A[r][tempCol] = Double.valueOf(s);
					tempCol+=1;
				}
				r+=1;
				c = tempCol;
			}
		} catch (NumberFormatException e) {
			// TODO Auto-generated catch block
			//e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		if (r > 3 || c >3)
		{
			System.out.println("Incorrect num of rows and columns(valid -less than 4)");
			return;
		}
		
		P1 obj=new P1();
		obj.getmatrix(A,r,c);
	}
	public void AmultiplicationNew(double[][] A,double[][] U, double[][] S, double[][] Vtranspose,int m,int n) {
		double Unew [][]= new double[m][n-1];
		double Snew [][]= new double[n-1][n-1];
		double Vtransposenew [][]= new double[n-1][n];
		
		int i,j = 0;
		
		//Unew
		for( i = 0; i < m; i++)
		{
			for( j = 0; j < n-1; j++)
			{
				Unew [i][j] =U [i][j];
				
			}
		}
		//Snew
		for( i = 0; i < n-1; i++)
		{
			for( j = 0; j < n-1; j++)
			{
				Snew [i][j] =S [i][j];
			}
		}
		//Vtransposenew
		for( i = 0; i < n-1; i++)
		{
			for( j = 0; j < n; j++)
			{
				Vtransposenew [i][j] =Vtranspose [i][j];
			}
		}
		System.out.println(" \nNew U after dim. reduction:");
		for( i = 0; i< m ; i++)
		{
			for( j = 0; j< n-1 ; j++)
				System.out.print(Math.round (Unew[i][j] * 10000.0) / 10000.0+ "\t");
					
			System.out.println();
		}
		System.out.println(" \nNew S after dim. reduction:");
		for( i = 0; i< n-1 ; i++)
		{
			for( j = 0; j< n-1 ; j++)
				System.out.print(Math.round (Snew[i][j] * 10000.0) / 10000.0+ "\t");
				
			System.out.println();
		}
		System.out.println("\nNew Vtranspose after dim. reduction:");
		for( i = 0; i< n-1 ; i++)
		{
			for( j = 0; j< n ; j++)
				System.out.print(Math.round (Vtransposenew[i][j] * 10000.0) / 10000.0+ "\t");
				
			System.out.println();
		}
				
		USVtransposeNewMult(A,Unew,Snew,Vtranspose,m,n);		
	}
	private void USVtransposeNewMult(double[][] A,double[][] unew, double[][] snew, double[][] vtranspose,int m,int n) {
		double UnewSnew [][]= new double[m][n-1];
		double USVtransposenew [][]= new double[m][n];
		int i,j;

		//Product of U and S matrices
		for( i = 0; i < m; i++)
		{
			for( j = 0; j < n-1; j++)
			{
				for(int k = 0; k < n-1; k++)
				{
					UnewSnew [i][j] += unew[i][k]*snew[k][j];
				}
			}
		}
		//Product of US and Vtranspose matrices
		for( i = 0; i < m; i++)
		{
			for( j = 0; j < n; j++)
			{
				for(int k = 0; k < n-1; k++)
				{
					USVtransposenew [i][j] += UnewSnew[i][k]*vtranspose[k][j];
				}
			}
		}
		System.out.println("\nProduct of matrix USVtranspose matrices after dim. reduction:");
		for( i = 0; i< m ; i++)
		{
			for( j = 0; j< n ; j++)
				System.out.print(Math.round (USVtransposenew[i][j] * 10000.0) / 10000.0+ "\t");
					
			System.out.println();
		}	
		Fnorm(A,unew,snew,vtranspose,USVtransposenew,m,n);
	}
	
	public void AmultiplicationNew1(double[][] A,double[][] U, double[][] S, double[][] Vtranspose,int m,int n) {
		double Unew [][]= new double[m][n-2];
		double Snew [][]= new double[n-2][n-2];
		double Vtransposenew [][]= new double[n-2][n];
		
		int i,j = 0;
		
		//Unew after 2nd dim
		for( i = 0; i < m; i++)
		{
			for( j = 0; j < n-2; j++)
			{
				Unew [i][j] =U [i][j];
				
			}
		}
		//Snew after 2nd dim
		for( i = 0; i < n-2; i++)
		{
			for( j = 0; j < n-2; j++)
			{
				Snew [i][j] =S [i][j];
			}
		}
		//Vtransposenew after 2nd dim
		for( i = 0; i < n-2; i++)
		{
			for( j = 0; j < n; j++)
			{
				Vtransposenew [i][j] =Vtranspose [i][j];
			}
		}
		System.out.println(" \nNew U after 2nd dim. reduction:");
		for( i = 0; i< m ; i++)
		{
			for( j = 0; j< n-2 ; j++)
				System.out.print(Math.round (Unew[i][j] * 10000.0) / 10000.0+ "\t");
					
			System.out.println();
		}
		System.out.println(" \nNew S after 2nd dim. reduction:");
		for( i = 0; i< n-2 ; i++)
		{
			for( j = 0; j< n-2 ; j++)
				System.out.print(Math.round (Snew[i][j] * 10000.0) / 10000.0+ "\t");
				
			System.out.println();
		}
		System.out.println("\nNew Vtranspose after 2nd dim. reduction:");
		for( i = 0; i< n-2 ; i++)
		{
			for( j = 0; j< n ; j++)
				System.out.print(Math.round (Vtransposenew[i][j] * 10000.0) / 10000.0+ "\t");
				
			System.out.println();
		}	
		USVtransposeNewMult1(A,Unew,Snew,Vtranspose,m,n);		
	}
	private void USVtransposeNewMult1(double[][] A,double[][] unew, double[][] snew, double[][] vtranspose,int m,int n) {
		double UnewSnew [][]= new double[m][n-2];
		double USVtransposenew [][]= new double[m][n];
		int i,j;

		//Product of U and S matrices
		for( i = 0; i < m; i++)
		{
			for( j = 0; j < n-2; j++)
			{
				for(int k = 0; k < n-2; k++)
				{
					UnewSnew [i][j] += unew[i][k]*snew[k][j];
				}
			}
		}
		//Product of US and Vtranspose matrices
		for( i = 0; i < m; i++)
		{
			for( j = 0; j < n; j++)
			{
				for(int k = 0; k < n-2; k++)
				{
					USVtransposenew [i][j] += UnewSnew[i][k]*vtranspose[k][j];
				}
			}
		}
		System.out.println("\nProduct of matrix USVtranspose matrices after 2nd dim. reduction:");
		for( i = 0; i< m ; i++)
		{
			for( j = 0; j< n ; j++)
				System.out.print(Math.round (USVtransposenew[i][j] * 10000.0) / 10000.0+ "\t");
					
			System.out.println();
		}	
		Fnormnew(A,unew,snew,vtranspose,USVtransposenew,m,n);
	}
	public void Fnormnew(double[][] a,double[][] unew, double[][] snew, double[][] vtranspose, double[][] uSVtranspose,int m,int n) {
		int i,j;
		double DiffAUSVtranspose [][]= new double[m][n];
		
		double ASqr = 0;
		for( i = 0; i < m; i++)
		{
			for( j = 0; j < n; j++)
			{
				DiffAUSVtranspose [i][j] = a[i][j]-uSVtranspose[i][j];				
			}
		}
		System.out.println("\nDifference of matrix original A and 2nd dim. reduction USVtranspose product matrix:");
		for( i = 0; i< m ; i++)
		{
			for( j = 0; j< n ; j++)
				System.out.print(Math.round (DiffAUSVtranspose[i][j] * 10000.0) / 10000.0+ "\t");
			
			System.out.println();
		}
		for( i = 0; i < m; i++)
		{
			for( j = 0; j < n; j++)
			{
				ASqr+=DiffAUSVtranspose[i][j] * DiffAUSVtranspose[i][j] ;						
			}
			
		}
		ASqr=Math.round (ASqr * 10000.0) / 10000.0 ;
		System.out.println("\nSquare sum of all terms :" +ASqr);
		double FN = Math.sqrt(ASqr);
		al.add(Math.round (FN * 10000.0) / 10000.0);
		if(m == 3 || n==3)
			{
				if (al.get(0) > al.get(1))
					System.out.println("\nBest Forbenius norm is after 2nd elimination :"+al.get(1));
				else
					System.out.println("\nBest Forbenius norm is after 1rst elimination :"+al.get(0));
			}
		else{
				if(al.size()==1)
					System.out.println("\nBest Forbenius norm is after 1rst elimination :"+al.get(0));
			}
	}
	public void Fnorm(double[][] a,double[][] unew, double[][] snew, double[][] vtranspose, double[][] uSVtranspose,int m,int n) {
		int i,j;
		double DiffAUSVtranspose [][]= new double[m][n];
		double ASqr = 0;
		for( i = 0; i < m; i++)
		{
			for( j = 0; j < n; j++)
			{
				DiffAUSVtranspose [i][j] = a[i][j]-uSVtranspose[i][j];
				DiffAUSVtranspose[i][j]= Math.round (DiffAUSVtranspose[i][j] * 10000.0) / 10000.0 ;
				
			}
		}
		System.out.println("\nDifference of matrix original A and dim. reduction USVtranspose product matrix:");
		for( i = 0; i< m ; i++)
		{
			for( j = 0; j< n ; j++)
				System.out.print(DiffAUSVtranspose[i][j]+ "\t");
			
			System.out.println();
		}
		for( i = 0; i < m; i++)
		{
			for( j = 0; j < n; j++)
			{
				ASqr+=DiffAUSVtranspose[i][j] * DiffAUSVtranspose[i][j] ;						
			}
			
		}
		ASqr=Math.round (ASqr * 10000.0) / 10000.0 ;
		System.out.println("\nSquare sum of all terms :" +ASqr);
		double FN = Math.sqrt(ASqr);
		al.add(Math.round (FN * 10000.0) / 10000.0);
			
		if((m == 3 || n==3) && (snew.length > 1))
			{	
				AmultiplicationNew1(a,unew,snew,vtranspose,m,n);
			}
		else{
				if(al.size()==1)
				System.out.println("\nBest Forbenius norm is after 1rst elimination :"+al.get(0));
			}
	}
	public void getmatrix(double A [][],int m,int n) {
		double AT [][]= new double[10][10];
		double ATA [][]= new double[10][10];
		int i,j = 0;
		
		System.out.println("Matrix A is: ");
		for( i = 0; i < m; i++)
		{
			for( j = 0; j < n; j++)
			{
				System.out.print(A[i][j]+"\t");
			}
			System.out.print("\n");
		}
		
		for( i = 0; i < n; i++)
		{
			for( j = 0; j < m; j++)
			{
				AT[i][j] = A[j][i];
			}
		}
		System.out.println("\nTranspose matrix AT of A is: ");
		for( i = 0; i < n; i++)
		{
			for( j = 0; j < m; j++)
			{
				System.out.print(AT[i][j]+"\t");
			}
			System.out.print("\n");
		}
		
		for( i = 0; i < n; i++)
		{
			for( j = 0; j < n; j++)
			{
				for(int k = 0; k < m; k++)
				{
					ATA [i][j] += AT[i][k]*A[k][j];
					ATA [i][j] =Math.round (ATA [i][j]  * 10000.0) / 10000.0 ;
				}
			}
		}
		System.out.println("\nProduct of A-transpose and A matrices:");
		for( i = 0; i< n ; i++)
		{
			for( j = 0; j< n ; j++)
				System.out.print(ATA[i][j]+ "\t");
			
			System.out.println();
		}
		
		Step2(A,ATA,n,m);
	}
	private void Step2(double[][] A,double[][] aTA,int n,int m) {
		double a = 0,c = 0, b = 0,d= 0, r1 = 0,r2 = 0,r3 = 0, rt1  = 0,rt2  = 0,rt3 = 0;
		int i=0, j=0;
		double S [][]= new double[n][n], Sinverse [][]= new double[n][n];
				
		 if (n==2)
	        {
	        	a = 1;
	        	b = -(aTA[0][0] + aTA[1][1]);
	        	c = -((aTA[0][1] * aTA[1][0])-(aTA[0][0] * aTA[1][1]));
	        	r3 = 0;
	        	
	        	System.out.println("\nCoefficients : a= " + a + " b= " + b + "c= " +c);
	        	double result = b * b - 4.0 * a * c;

	            if (result > 0.0) {
	                 r1 = (-b + Math.pow(result, 0.5)) / (2.0 * a);
	                 r2 = (-b - Math.pow(result, 0.5)) / (2.0 * a);
	                System.out.println("\nThe roots are " + Math.round (r1 * 10000.0) / 10000.0 + " and " + Math.round (r2 * 10000.0) / 10000.0);
	            } else if (result == 0.0) {
	                 r1 = -b / (2.0 * a);
	                System.out.println("\nThe root is " + Math.round (r1 * 10000.0) / 10000.0);
	            } else {
	                System.out.println("\nThe equation has no real roots.");
	            }
	        }
		if (n==2)
		   {
		        S[0][0] = Math.sqrt(r1);
		        S[0][0] = Math.round (S[0][0]  * 10000.0) / 10000.0 ;
		        S[0][1] = 0;
		        S[1][0] = 0;
		        S[1][1] = Math.sqrt(r2);
		        S[1][1] = Math.round (S[1][1]  * 10000.0) / 10000.0 ;
		        System.out.println("\nMatrix S is(Singular values):");
		  for( i = 0; i< n ; i++)
		  {
			for( j = 0; j< n ; j++)
					System.out.print(S[i][j]+ "\t");
				
				System.out.println();
		  }
	      	 for( i = 0; i< n ; i++)
			for( j = 0; j< n ;j++)
			{
				if (i!=j)
				{
					Sinverse[i][j]=0;}
				else
				{
					Sinverse [i][j]= 1 / S[i][j];
				}
			}
	      
		System.out.println("\nInverse of S is:");
		for( i = 0; i< n ; i++)
		{
			for( j = 0; j< n ; j++)
				System.out.print(Math.round (Sinverse[i][j] * 10000.0) / 10000.0+ "\t");
				
			System.out.println();
		}
		Step4V(A,aTA,S,Sinverse,r1,r2,r3,n,m);
	       }
	     if(n==3)
	        {
	        	a= -1;
	        	b= aTA[0][0] + -(-aTA[2][2] + -aTA[1][1]);
	    		c= ( -(aTA[1][1] * aTA[2][2] - aTA[2][1]*aTA[1][2])) + aTA[0][0]*(-aTA[1][1]+ -aTA[2][2]) +
	    				(-aTA[0][1] *-aTA[1][0])+ aTA[0][2] * aTA[2][0] ;
	    		d= aTA[0][0]*(aTA[1][1] * aTA[2][2] - aTA[2][1] * aTA[1][2]) + -aTA[0][1]*
	    				((aTA[1][0] * aTA[2][2]) - (aTA[1][2] * aTA[2][0]))+
	    				aTA[0][2]*(aTA[1][0] * aTA[2][1] - aTA[1][1] * aTA[2][0]);
	    		System.out.println("\nCoefficients : a= " + a + " b= " + b + " c= " +c + " d= " +d);
	    		{
	    		if (a == 0.0)
	    			{
	    				throw new RuntimeException ("Cubic.solve(): a = 0");
	    			}
	    		double denom = a;a = b/denom;b = c/denom;c = d/denom;
	    		int nRoots;
	    			 
	    		// Commence solution.
	    		double a_over_3 = a / 3.0;double Q = (3*b - a*a) / 9.0;
	    		double Q_CUBE = Q*Q*Q;double R = (9*a*b - 27*c - 2*a*a*a) / 54.0;
	    		double R_SQR = R*R;double D = Q_CUBE + R_SQR;

	    		if (D < 0.0)
	    			{
	    			// Three unequal real roots.
	    			nRoots = 3;
	    			double theta = Math.acos (R / Math.sqrt (-Q_CUBE));
	    			double SQRT_Q = Math.sqrt (-Q);
	    			rt1 = 2.0 * SQRT_Q * Math.cos (theta/3.0) - a_over_3;
	    			rt1= Math.round (rt1 ) ;
	    			rt2 = 2.0 * SQRT_Q * Math.cos ((theta+TWO_PI)/3.0) - a_over_3;
	    			rt2= Math.round (rt2 );
	    			rt3 = 2.0 * SQRT_Q * Math.cos ((theta+FOUR_PI)/3.0) - a_over_3;
	    			rt3= Math.round (rt3 );
	    			}
	    		else if (D > 0.0)
	    			{
	    			// One real root.
	    			nRoots = 1;
	    			double SQRT_D = Math.sqrt (D);
	    			double S1 = Math.cbrt (R + SQRT_D);
	    			double T = Math.cbrt (R - SQRT_D);
	    			rt1 = (S1 + T) - a_over_3;
	    			rt1= Math.round (rt1 * 10000.0) / 10000.0 ;
	    			rt2 = Double.NaN;
	    			rt3 = Double.NaN;
	    			}
	    		else
	    			{
	    			// Three real roots, at least two equal.
	    			nRoots = 3;
	    			double CBRT_R = Math.cbrt (R);
	    			rt1 = 2*CBRT_R - a_over_3;
	    			rt1= Math.round (rt1 * 10000.0) / 10000.0 ;
	    			rt2 = rt3 = CBRT_R - a_over_3;
	    			rt2= Math.round (rt2 * 10000.0) / 10000.0 ;
	    			}
	    		if (rt1 < rt2){double tmp = rt1; rt1 = rt2; rt2 = tmp;}
	    		if (rt2 < rt3){double tmp = rt2; rt2 = rt3; rt3 = tmp;}
	    		if (rt1 < rt2){double tmp = rt1; rt1 = rt2; rt2 = tmp;}
	    		System.out.println("\nroot 1 = " + rt1);
	    		if (nRoots == 3) {
	    			System.out.println("root 2 = " + rt2);
	    			System.out.println("root 3 = " + rt3);}
	    		}
	        }
	        
	        if(n==3) {
	        	S[0][0] = Math.sqrt(rt1);
		        S[0][0] = Math.round (S[0][0]  * 10000.0) / 10000.0 ;
		        S[0][1] = 0;
		        S[1][0] = 0;
		        S[1][1] = Math.sqrt(rt2);
		        S[1][1] = Math.round (S[1][1]  * 10000.0) / 10000.0 ;
		        S[0][2] = 0;
		        S[2][0] = 0;
		        S[2][1] = 0;
		        S[2][2] = Math.sqrt(rt3);
		        S[2][2] = Math.round (S[2][2]  * 10000.0) / 10000.0 ;
		        System.out.println("\nMatrix S is(Singular values):");
		        for( i = 0; i< n ; i++)
				{
				for( j = 0; j< n ; j++)
					System.out.print(S[i][j]+ "\t");
				
				System.out.println();
				}
		        for( i = 0; i< n ; i++)
					for( j = 0; j< n ;j++)
					{
						if (i!=j)
						{
							Sinverse[i][j]=0;}
						else
						{
						Sinverse [i][j]=1 / S[i][j];
						Sinverse[i][j]= Math.round (Sinverse[i][j] * 10000.0) / 10000.0 ;
						}
					}
				System.out.println("\nInverse of S is :");
				for( i = 0; i< n ; i++)
					{
					for( j = 0; j< n ; j++)
						System.out.print(Sinverse[i][j]+ "\t");
					
					System.out.println();
					}
				
			Step4V(A,aTA,S,Sinverse,rt1,rt2,rt3,n,m);
			}
	}
	private void Step4V(double[][] A,double[][] aTA,double[][] s, double[][] sinverse, double r1,double r2,double r3, int n,int m) {
		
	double aTAnew [][]= new double[n][n];
	double aTAnew2 [][]= new double[n][n];
	double aTAnew3 [][]= new double[n][n];
	double V [][] = new double[n][n];
	double VT [][]= new double[n][n];
	double X1 [][]= new double[n][1];
	double X2 [][]= new double[n][1];
	double X3 [][]= new double[n][1];
	int i,j;
		
	aTAnew[0][0] = aTA[0][0] - r1;aTAnew[0][1] = aTA[0][1];
        aTAnew[1][0] = aTA[1][0];aTAnew[1][1] = aTA[1][1] - r1;
        
        aTAnew2[0][0] = aTA[0][0] - r2; aTAnew2[0][1] = aTA[0][1];       
        aTAnew2[1][0] = aTA[1][0];aTAnew2[1][1] = aTA[1][1] - r2;
        
        if(n==3) {
        	aTAnew3[0][0] = aTA[0][0] - r3;aTAnew3[0][1] = aTA[0][1];
        	aTAnew3[1][0] = aTA[1][0];aTAnew3[1][1] = aTA[1][1] - r3;
        }
        
        if(n==3){
        	aTAnew[2][2] = aTA[2][2] - r1;
            aTAnew[0][2] = aTA[0][2];aTAnew[1][2] = aTA[1][2];
            aTAnew[2][1] = aTA[2][1]; aTAnew[2][0] = aTA[2][0];
            
            aTAnew2[2][2] = aTA[2][2] - r2;
            aTAnew2[0][2] = aTA[0][2];aTAnew2[1][2] = aTA[1][2];
            aTAnew2[2][1] = aTA[2][1];aTAnew2[2][0] = aTA[2][0];
            
            aTAnew3[2][2] = aTA[2][2] - r3; aTAnew3[0][2] = aTA[0][2];
            aTAnew3[1][2] = aTA[1][2];aTAnew3[2][1] = aTA[2][1];
            aTAnew3[2][0] = aTA[2][0];
        }
        
        if(n==2) {
        	double e,x,y,r;
        	e = r1;
        	x = -aTA[0][1]; y = (aTA[0][0]-e);
        	r = Math.sqrt(x*x+y*y);
	    
	    if( r > 0) { 
	    	x =x/r; y= y/ r;}
	    else {
	        x = e-aTA[1][1]; y = aTA[1][0];
	        r = Math.sqrt(x*x+y*y);
	        if( r > 0) { x /= r; y /= r;}
	        else {
	            x = 1; y = 0;
	           }
	    }
	    
	    X1[0][0]=x;
	    X1[1][0]=y;
	    
	    	e = r2;
	    	x = -aTA[0][1]; y = (aTA[0][0]-e);
	    	r = Math.sqrt(x*x+y*y);
	    if( r > 0) {
	    	x =x/r; y = y/r; 
	    }
	    else {
	        x = e-aTA[1][1]; y = aTA[1][0];
	        r = Math.sqrt(x*x+y*y);
	        if( r > 0) { x /= r; y /= r; }
	        else {
	            x = 0; y = 1;
	        }
	    }
	    
	    X2[0][0]=x;
	    X2[1][0]=y;
        }
		        
        if(n==3) {
        	double e=-aTAnew[0][0]; double a= aTAnew[0][1];
        	double b= aTAnew[0][2]; double f=-aTAnew[1][0];
        	double c= aTAnew[1][1]; double d= aTAnew[1][2];
        	int x=1;
        	double det = ((a) * (d) - (b) * (c));
        	double y = ((d) * (e) - (b) * (f)) / det;
            double z = ((a) * (f) - (c) * (e)) / det;
    
            X1[0][0]=Math.round(x)/Math.sqrt(x*x+y*y+z*z);
    	    X1[1][0]=Math.round(y)/Math.sqrt(x*x+y*y+z*z);
    	    X1[2][0]=Math.round(z)/Math.sqrt(x*x+y*y+z*z);
    	    
    	    double e1=-aTAnew2[0][0];double a1= aTAnew2[0][1];
        	double b1= aTAnew2[0][2];double f1=-aTAnew2[1][0];
        	double c1= aTAnew2[1][1];double d1= aTAnew2[1][2];
        	int x1=1;
        	double det1 = ((a1) * (d1) - (b1) * (c1));  
            double y1 = ((d1) * (e1) - (b1) * (f1)) / det1;
            double z1 = ((a1) * (f1) - (c1) * (e1)) / det1;
    
            X2[0][0]=Math.round(x1)/Math.sqrt(x1*x1+y1*y1+z1*z1);
    	    X2[1][0]=Math.round(y1)/Math.sqrt(x1*x1+y1*y1+z1*z1);
    	    X2[2][0]=Math.round(z1)/Math.sqrt(x1*x1+y1*y1+z1*z1);
    	    
    	    double e2=-aTAnew3[0][0]; double a2= aTAnew3[0][1];
        	double b2= aTAnew3[0][2]; double f2=-aTAnew3[1][0];
        	double c2= aTAnew3[1][1]; double d2= aTAnew3[1][2];
        	int x2=1;
        	double det2 = ((a2) * (d2) - (b2) * (c2)); 
            double y2 = ((d2) * (e2) - (b2) * (f2)) / det2;
            double z2 = ((a2) * (f2) - (c2) * (e2)) / det2;
    
            X3[0][0]=Math.round(x2)/Math.sqrt(x2*x2+y2*y2+z2*z2);
    	    X3[1][0]=Math.round(y2)/Math.sqrt(x2*x2+y2*y2+z2*z2);
    	    X3[2][0]=Math.round(z2)/Math.sqrt(x2*x2+y2*y2+z2*z2);
        }
        
        for( i = 0; i < n; i++)
		{
			for( j = 0; j < 1; j++)
			{
				V[i][j] = X1[i][j];
			}
		}
        
		for( i = 0; i < n; i++)
		{
			for( j = 0; j < 1; j++)
			{
				V[i][j+1] = X2[i][j];
			}
		}
		if(n==3) {
		for( i = 0; i < n; i++)
		{
			for( j = 0; j < 1; j++)
			{
				V[i][j+2] = X3[i][j];
			}
		} }
		System.out.println("\nMatrix V is(Right singular vectors): ");
		for( i = 0; i < n; i++)
		{
			for( j = 0; j < n; j++)
			{
				System.out.print(Math.round (V[i][j] * 10000.0) / 10000.0 +"\t");
			}
			System.out.print("\n");
		}
		
		for( i = 0; i < n; i++)
		{
			for( j = 0; j < n; j++)
			{
				VT[i][j] = V[j][i];
			}
		}
		System.out.println("\nTranspose matrix VT of V is: ");
		for( i = 0; i < n; i++)
		{
			for( j = 0; j < n; j++)
			{
				System.out.print(Math.round (VT[i][j] * 10000.0) / 10000.0+"\t");
			}
			System.out.print("\n");
		}
		Umultplication(A,V,VT,s,sinverse,n,m);
			
	}
	private void Umultplication(double[][] A, double[][] V,double[][] VT,double[][] s, double[][] sinverse,int n,int m) {
		double AV [][]= new double[m][n];
		double U [][]= new double[m][n];
		
		int i,j = 0;
		//Product of A and V matrices
		for( i = 0; i < m; i++)
		{
			for( j = 0; j < n; j++)
			{
				for(int k = 0; k < n; k++)
				{
					//AAT [i][j] += AT[k][i]*A[j][k];
					AV [i][j] +=A[i][k]*V[k][j];
				}
			}
		}
		//Product of AV and Sinverse matrices
		for( i = 0; i < m; i++)
		{
			for( j = 0; j < n; j++)
			{
				for(int k = 0; k < n; k++)
				{
					//AAT [i][j] += AT[k][i]*A[j][k];
					U [i][j] += AV[i][k]*sinverse[k][j];
				}
			}
		}
		System.out.println("\nMatrix U(Left singular vectors):");
		for( i = 0; i< m ; i++)
			{
			for( j = 0; j< n ; j++)
				System.out.print(Math.round (U[i][j] * 10000.0) / 10000.0+ "\t");
			
			System.out.println();
			}
		Amultiplication(A,U,s,VT,m,n);
	}
	
	public void Amultiplication(double[][] A,double[][] U, double[][] S, double[][] Vtranspose,int m,int n) {
		double US [][]= new double[m][n];
		double USVtranspose [][]= new double[m][n];
		
		int i,j = 0;
		//Product of U and S matrices
		for( i = 0; i < m; i++)
		{
			for( j = 0; j < n; j++)
			{
				for(int k = 0; k < n; k++)
				{
					US [i][j] += U[i][k]*S[k][j];
				}
			}
		}
		System.out.println("\nMatrix product of U and S:");
		for( i = 0; i< m ; i++)
		{
			for( j = 0; j< n ; j++)
				System.out.print(Math.round (US[i][j] * 10000.0) / 10000.0+ "\t");
			
			System.out.println();
		}
		//Product of US and Vtranspose matrices
		for( i = 0; i < m; i++)
		{
			for( j = 0; j < n; j++)
			{
				for(int k = 0; k < n; k++)
				{
					USVtranspose [i][j] += US[i][k]*Vtranspose[k][j];
					USVtranspose[i][j]= Math.round (USVtranspose[i][j] * 10000.0) / 10000.0 ;
				}
			}
		}
		System.out.println("\nMatrix A = U*S*V^T:");
		for( i = 0; i< m ; i++)
		{
			for( j = 0; j< n ; j++)
				System.out.print(Math.round (USVtranspose[i][j] * 10000.0) / 10000.0+ "\t");
			
			System.out.println();
		}
		AmultiplicationNew(A,U,S,Vtranspose,m,n);
	}
}
