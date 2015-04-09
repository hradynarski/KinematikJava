package kinematik;

public class LA {  //linear algebra

	public static double dist(double[] a, double[] b) {
		double r = 0;
		for (int i = 0; i < a.length; i++)
			r += (a[i] - b[i]) * (a[i] - b[i]);
		return Math.sqrt(r);
	}

	public static double[] sub(double[] a, double[] b) {
		double[] re = new double[a.length];
		for (int i = 0; i < a.length; i++)
			re[i] = a[i] - b[i];
		return re;
	}

	public static double mag(double[] a) {
		double r = 0;
		for (int i = 0; i < a.length; i++)
			r += a[i] * a[i];
		r=Math.sqrt(r);
		return r;
	}
	
	public static double[] normalize(double[] a) {
		double r = mag(a);
		double[] re=new double[a.length];
		for (int i = 0; i < a.length; i++)
			re[i]= a[i]/r;
		return re;
	}
	
	public static double dot(double[] a, double[] b) {
		double re = 0;
		for (int i = 0; i < a.length; i++)
			re+= a[i] * b[i];
		return re;
	}
	
	public static double[] mul(double[] a, double s) {
		double[] re=new double[a.length];
		for (int i = 0; i < a.length; i++)
			re[i]= a[i] *s;
		return re;
	}
	
	public static double[] cross(double[] b, double[] c) {
		double[] re = new double[b.length];
		re[0] = b[1] * c[2] - b[2] * c[1];
		re[1] = b[2] * c[0] - b[0] * c[2];
		re[2] = b[0] * c[1] - b[1] * c[0];
		return re;
	}
}
