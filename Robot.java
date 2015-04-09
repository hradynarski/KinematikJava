package kinematik;


public class Robot {
	private static final double PI = Math.PI;
	private static final double HP = 0.5 * Math.PI;
	final double[] a,d;
	private final double l2, ad2;
	
	public Robot(double[] a, double[] d) { // a,d: DH parameters, corresponding to the geometry of your robot
		this.a=a;
		this.d=d;
		l2 = Math.sqrt(a[2] * a[2] + d[3] * d[3]);
		ad2 = Math.atan2(a[2], d[3]);
	}
	
	public double[] inverse(double[] v, Double A4_deg, double[][] base, double[][] tool) {
		double[][] Pos = matrix(v[0], v[1], v[2], v[3], v[4], v[5]);
		double[][] goal;
		if(base==null)
			goal=Pos;
		else
			goal= mul34(base, Pos); // in WORLD frame BASE*POS(in BASE)= GOAL
		double[][] T6;
		if (tool == null) {
			T6 = goal;
		} else {
			double[][] intool = inverse34(tool);
			T6 = mul34(goal, intool); // T6*TOOL = GOAL= BASE*POS(in BASE)
		}
		double[] inreDeg = inverse(T6, A4_deg);// (T6, thetaDeg[3]);
		return inreDeg;
	}
	
	public double[][] forward(double[] degs) {  
		double[] ts= deg_rad(degs);
		double[] c=new double[6];
		double[] s=new double[6];
		for(int i=0;i<6;i++){
			c[i]= Math.cos(ts[i]);
			s[i]= Math.sin(ts[i]);
		}
		double[][] m123=new double[3][];
		m123[0] = new double[] { c[0] * (c[1] * c[2] - s[1] * s[2]),     s[0],   c[0] * (c[1] * s[2] + s[1] * c[2]),      c[0] * (a[2] * (c[1] * c[2] - s[1] * s[2]) + a[1] * c[1]) + a[0] * c[0] };
		m123[1] = new double[] { s[0]*  (c[1] * c[2] - s[1] * s[2]),   -c[0],    s[0] * (c[1] * s[2] + s[1] * c[2]),      s[0]* (a[2] * (c[1] * c[2] - s[1] * s[2]) + a[1] * c[1]) + a[0] * s[0]};
		m123[2] = new double[] { s[1] * c[2] +c[1] * s[2],                0,    s[1] * s[2] - c[1] * c[2],                 a[2] * (s[1] * c[2] + c[1] * s[2]) + a[1] * s[1] +d[0]};
		double[][] m456=new double[3][];
		m456[0] = new double[] {c[3]*c[4]*c[5]-s[3]*s[5],   -c[3]*c[4]*s[5]-s[3]*c[5]     , c[3] * s[4],     c[3]*s[4]*d[5] };
		m456[1] = new double[] {s[3]*c[4]*c[5]+c[3]*s[5],   -s[3]*c[4]*s[5]+c[3]*c[5]     , s[3] * s[4],     s[3]*s[4]*d[5] };
		m456[2] = new double[] { -s[4]*c[5],                  s[4]*s[5],                     c[4],             c[4]*d[5]+d[3]};
		double[][] arr = mul34(m123, m456);
		return arr;
	}
	
	public double[] forward(double[]  degs, double[][] base, double[][] tool) {
		double[][] T6 = forward(degs);
		double[][] goal;
		if (tool == null) {
			goal = T6;
		} else {
			goal = mul34(T6, tool);
		}
		double[][] pos;
		if (base == null) {
			pos = goal;
		} else {
			double[][] inbase = inverse34(base);
			pos = mul34(inbase, goal);
		}
		double[] as = Robot.ABC(pos);
		return new double[] { pos[0][3], pos[1][3], pos[2][3], as[0], as[1], as[2] };
	}

	public double[] inverse(double[][] T6, Double A4_deg) {
		double[] theta=new double[9]; 
		double[] center= mul34(T6, new double[]{0,0, -d[5]});
		theta[0] = Math.atan2(center[1], center[0]); // or -atan2     choice one possibility

		double ll = Math.sqrt(center[0] * center[0] + center[1] * center[1]);
		double[] p1 = { a[0] * center[0] / ll, a[0] * center[1] / ll, d[0] };
		double l3 =LA.dist(center, p1);
		double l1 = a[1];
		double beta = Math.acos((l1 * l1 + l3 * l3 - l2 * l2) / (2 * l1 * l3));
		double ttl = Math.sqrt((center[0] - p1[0]) * (center[0] - p1[0]) + (center[1] - p1[1]) * (center[1] - p1[1]));
		if (p1[0] * (center[0] - p1[0]) < 0) // opposite side
			ttl = -ttl;
		double al = Math.atan2(center[2] - p1[2], ttl);
		theta[1] =beta+al; // choice one possibility
		double gama = Math.acos((l1 * l1 + l2 * l2 - l3 * l3) / (2 * l1 * l2));
		theta[2] = gama - ad2 - HP;

		double[][] arr = new double[4][];
		double[] c = new double[3];
		double[] s=new double[3];
		for(int i=0;i<3;i++){
			c[i]= Math.cos(theta[i]);
			s[i]= Math.sin(theta[i]);
		}
		arr[0] = new double[] { c[0] * (c[1] * c[2] - s[1] * s[2]),     s[0],   c[0] * (c[1] * s[2] + s[1] * c[2]),      c[0] * (a[2] * (c[1] * c[2] - s[1] * s[2]) + a[1] * c[1]) + a[0] * c[0] };
		arr[1] = new double[] { s[0]*  (c[1] * c[2] - s[1] * s[2]),   -c[0],    s[0] * (c[1] * s[2] + s[1] * c[2]),      s[0]* (a[2] * (c[1] * c[2] - s[1] * s[2]) + a[1] * c[1]) + a[0] * s[0]};
		arr[2] = new double[] { s[1] * c[2] +c[1] * s[2],                0,    s[1] * s[2] - c[1] * c[2],                 a[2] * (s[1] * c[2] + c[1] * s[2]) + a[1] * s[1] +d[0]};
		double[][] in123= inverse34(arr);
		double[][] mr = mul34(in123, T6);
		double c5 = mr[2][2];
		if (Math.abs(c5 - 1) < 0.000001) { //singularity
			double A4=-PI*A4_deg/180; //see deg_rad
			double c4 = Math.cos(A4);
			double s4 = Math.sin(A4);
			double s6 = c4 * mr[1][0] - s4 * mr[0][0];
			double c6;
			if (Math.abs(c4) > Math.abs(s4))
				c6 = (mr[0][0] + s4 * s6) / c4;
			else
				c6 = (mr[1][0] - c4 * s6) / s4;
			theta[3] = A4;
			theta[4] = 0;
			theta[5] =   Math.atan2(s6, c6);
			if (Math.abs(c6) > 1 || Math.abs(s6) > 1)
				throw new RuntimeException();
		} else {
			double ang = Math.atan2(mr[1][2], mr[0][2]);
			theta[3] = ang;
			theta[4] = Math.acos(c5); // *********
			theta[5] = Math.atan2(mr[2][1], -mr[2][0]);
		}
		double[] inreDeg = rad_deg(theta);
		return inreDeg;
	}

	private double[] deg_rad(double[] ds) {
		double[] rd = new double[6];
		for (int i = 0; i < 6; i++)
			rd[i] = ds[i] * PI / 180;
		rd[2] -= HP;
		rd[5] += PI;
		for (int i = 0; i < 6; i++)
			rd[i] = -rd[i];
		return rd;
	}

	private double[] rad_deg(double[] ds) {
		double[] rd = new double[6];
		for (int i = 0; i < 6; i++)
			rd[i] = -ds[i];
		rd[2] += HP;
		rd[5] -= PI;
		for (int i = 0; i < 6; i++)
			rd[i] = rd[i] * 180 / PI;
		return rd;
	}
	
	public static double[][] mul34(double[][] a, double[][] b) { // multiple two 3*4 matrices
		double[][] re = new double[3][4];
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 4; j++) {
				double b3j = (j == 3 ? 1 : 0);
				re[i][j] = a[i][0] * b[0][j] + a[i][1] * b[1][j] + a[i][2] * b[2][j] + a[i][3] * b3j;
			}
		}
		return re;
	}

	public static double[] mul34(double[][] a, double[] b) { // multiple a 3*4 matrix with a vector
		double[] re = new double[3];
		for (int i = 0; i < 3; i++)
			re[i] = a[i][0] * b[0] + a[i][1] * b[1] + a[i][2] * b[2] + a[i][3];
		return re;
	}
	
	public static double[][] matrix(double x, double y, double z, double aDeg, double bDeg, double cDeg) {
		double a = -aDeg * PI / 180;
		double b = -bDeg * PI / 180;
		double c = -cDeg * PI / 180;
		double ca = Math.cos(a);
		double sa = Math.sin(a);
		double cb = Math.cos(b);
		double sb = Math.sin(b);
		double cc = Math.cos(c);
		double sc = Math.sin(c);
		double[][] tt = new double[3][];
		tt[0] = new double[] { ca * cb, sa * cc + ca * sb * sc, sa * sc - ca * sb * cc, x };
		tt[1] = new double[] { -sa * cb, ca * cc - sa * sb * sc, ca * sc + sa * sb * cc, y };
		tt[2] = new double[] { sb, -cb * sc, cb * cc, z };
		return tt;
	}
	
	public static double[][] matrix(double aDeg, double bDeg, double cDeg) {
		double a = -aDeg * PI / 180;
		double b = -bDeg * PI / 180;
		double c = -cDeg * PI / 180;
		double ca = Math.cos(a);
		double sa = Math.sin(a);
		double cb = Math.cos(b);
		double sb = Math.sin(b);
		double cc = Math.cos(c);
		double sc = Math.sin(c);
		double[][] tt = new double[3][];
		tt[0] = new double[] { ca * cb, sa * cc + ca * sb * sc, sa * sc - ca * sb * cc };
		tt[1] = new double[] { -sa * cb, ca * cc - sa * sb * sc, ca * sc + sa * sb * cc};
		tt[2] = new double[] { sb, -cb * sc, cb * cc };
		return tt;
	}
	
	public static double[] PosByVect(double[] _dx, double[] _dxy) {//_dx:  x-axis    _dxy: a vector  on XY plane,     
		double[] dx= LA.normalize(_dx);
		double[] tt= LA.mul(dx, LA.dot(_dxy, dx));
		double[] _dy= LA.sub(_dxy, tt);
		double[] dy= LA.normalize(_dy);
		double[] dz= LA.cross(dx, dy);
				
		double cacb=dx[0];
		double sacb=-dx[1];
		double sb=dx[2];
		double cbsc=-dy[2];
		double cbcc= dz[2];
	
		double cb = Math.sqrt(1 - sb * sb); // + -  similar to ABC()
		double a = Math.atan2(sacb, cacb) * -180 /PI;
		double b = Math.atan2(sb, cb) * -180 / PI;
		double c = Math.atan2(cbsc, cbcc) * -180 / PI;
		return new double[] { a, b, c };
	}

	public static double[] ABC(double[][] m) { //Euler angles from 3*3 matrix
		double sb = m[2][0];
		double cb = Math.sqrt(1 - sb * sb); // + -
		double ca = m[0][0];
		double sa = -m[1][0];
		double cc = m[2][2];
		double sc = -m[2][1];
		double a = Math.atan2(sa, ca) * -180 /PI;
		double b = Math.atan2(sb, cb) * -180 / PI;
		double c = Math.atan2(sc, cc) * -180 / PI;
		return new double[] { a, b, c };
	}
	public static double[] flipABC( double[] abc){
		return flipABC(abc[0], abc[1], abc[2]);
	}
	public static double[] flipABC( double a, double b, double c){
		double na= a>0? (a-180):(a+180);
		double nb= b>0?  (180-b):(-180-b);
		double nc= c>0? (c-180):(c+180);
		return new double[]{na, nb, nc};
	}
	public static double[][] inverse34(double[][] m) { // det=1,  row 3 is {0,0,0,1}
		double[][] v = new double[3][4];
		v[0][0] = -m[1][2] * m[2][1] + m[1][1] * m[2][2];
		v[0][1] = m[0][2] * m[2][1] - m[0][1] * m[2][2];
		v[0][2] = -m[0][2] * m[1][1] + m[0][1] * m[1][2];
		v[0][3] = m[0][3] * m[1][2] * m[2][1] - m[0][2] * m[1][3] * m[2][1] - m[0][3] * m[1][1] * m[2][2] + m[0][1] * m[1][3] * m[2][2] + m[0][2] * m[1][1] * m[2][3] - m[0][1]
				* m[1][2] * m[2][3];
		v[1][0] = m[1][2] * m[2][0] - m[1][0] * m[2][2];
		v[1][1] = -m[0][2] * m[2][0] + m[0][0] * m[2][2];
		v[1][2] = m[0][2] * m[1][0] - m[0][0] * m[1][2];
		v[1][3] = m[0][2] * m[1][3] * m[2][0] - m[0][3] * m[1][2] * m[2][0] + m[0][3] * m[1][0] * m[2][2] - m[0][0] * m[1][3] * m[2][2] - m[0][2] * m[1][0] * m[2][3] + m[0][0]
				* m[1][2] * m[2][3];
		v[2][0] = -m[1][1] * m[2][0] + m[1][0] * m[2][1];
		v[2][1] = m[0][1] * m[2][0] - m[0][0] * m[2][1];
		v[2][2] = -m[0][1] * m[1][0] + m[0][0] * m[1][1];
		v[2][3] = m[0][3] * m[1][1] * m[2][0] - m[0][1] * m[1][3] * m[2][0] - m[0][3] * m[1][0] * m[2][1] + m[0][0] * m[1][3] * m[2][1] + m[0][1] * m[1][0] * m[2][3] - m[0][0]
				* m[1][1] * m[2][3];
		return v;
	}

}