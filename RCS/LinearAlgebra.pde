public class Matrix {
  public int n, m;
  public double[][] entries;
  
  public Matrix() {
    n = 0;
    m = 0;
  }
  
  public Matrix(int n, int m) {
    this.n = n;
    this.m = m;
    entries = new double[n][];
    for(int i = 0; i < n; i++) {
      entries[i] = new double[m];
    }
  }
  
  public Matrix(Matrix m) {
    this(m.entries);
    this.n = m.n;
    this.m = m.m;
  }
  
  public Matrix(double[][] entriesToCopy) {
    n = entriesToCopy.length;
    m = entriesToCopy[0].length;
    entries = new double[entriesToCopy.length][];
    for(int i = 0; i < entries.length; i++) {
      entries[i] = new double[entriesToCopy[i].length];
      for(int j = 0; j < entries[i].length; j++) {
        entries[i][j] = entriesToCopy[i][j];
      }
    }
  }
  
  public Matrix(int n, int m, double... values) {
    int index = 0;
    this.n = n;
    this.m = m;
    entries = new double[n][];
    for(int i = 0; i < n; i++) {
      entries[i] = new double[m];
      for(int j = 0; j < m; j++) {
        entries[i][j] = values[index];
        index++;
      }
    }
  }
  
  public Matrix Transpose() {
    Matrix m = new Matrix(this.m, this.n);
    double[][] result = m.entries;
    for(int i = 0; i < result.length; i++) {
      result[i] = new double[entries.length];
      for(int j = 0; j < result[i].length; j++) {
        result[i][j] = entries[j][i];
      }
    }
    return m;
  }
  
  public double Determinant() {
    if(n == 1 && m == 1) return entries[0][0];
    if(n == 2 && m == 2) {
      return entries[0][0] * entries[1][1] - entries[0][1] * entries[1][0];
    }
    if(n == 3 && m == 3) {
      double a = entries[0][0]; double b = entries[0][1]; double c = entries[0][2];
      double d = entries[1][0]; double e = entries[1][1]; double f = entries[1][2];
      double g = entries[2][0]; double h = entries[2][1]; double i = entries[2][2];
      return a*(e*i-f*h)-b*(d*i-f*g)+c*(d*h-e*g);
    }
    return 0;
  }
  
  public Matrix Inverse() {
    if(n != m) return null;
    if(n == 0) return new Matrix();
    if(n == 1) return new Matrix(1,1,1 / entries[0][0]);
    if(n == 2) {
      double invDeterminant = 1 / Determinant();
      return new Matrix(2,2,
         entries[1][1] * invDeterminant, -entries[0][1] * invDeterminant,
        -entries[1][0] * invDeterminant,  entries[0][0] * invDeterminant);
    }
    if(n == 3) {
      double a = entries[0][0]; double b = entries[0][1]; double c = entries[0][2];
      double d = entries[1][0]; double e = entries[1][1]; double f = entries[1][2];
      double g = entries[2][0]; double h = entries[2][1]; double i = entries[2][2];
      double invD = 1 / Determinant();
      
      return new Matrix(3,3,
        invD*(e*i-f*h), invD*(c*h-b*i), invD*(b*f-c*e),
        invD*(f*g-d*i), invD*(c*g-a*i), invD*(c*d-a*f),
        invD*(d*h-e*g), invD*(b*g-a*h), invD*(a*e-b*d));
    }
    
    return null;
  }
  
  public Matrix add(Matrix m2) {
    Matrix result = new Matrix(this.n, this.m);
    for(int i = 0; i < n; i++) {
      for(int j = 0; j < m; j++) {
        result.entries[i][j] = this.entries[i][j] + m2.entries[i][j];
      }
    }
    return result;
  }
  
  public Matrix sub(Matrix m2) {
    Matrix result = new Matrix(this.n, this.m);
    for(int i = 0; i < n; i++) {
      for(int j = 0; j < m; j++) {
        result.entries[i][j] = this.entries[i][j] - m2.entries[i][j];
      }
    }
    return result;
  }
  
  public Matrix mult(Matrix right) {
    Matrix result = new Matrix(this.n, right.m);
    for(int i = 0; i < result.entries.length; i++) {
      for(int j = 0; j < result.entries[i].length; j++) {
        double sum = 0;
        for(int k = 0; k < this.m; k++) {
          sum += entries[i][k] * right.entries[k][j];
        }
        result.entries[i][j] = sum;
      }
    }
    return result;
  }
  
  public Matrix mult(double scalar) {
    Matrix result = new Matrix(n, m);
    for(int i = 0; i < n; i++) {
      for(int j = 0; j < m; j++) {
        result.entries[i][j] = entries[i][j] * scalar;
      }
    }
    return result;
  }
  
  public void printMatrix() {
    for(int i = 0; i < n; i++) {
      for(int j = 0; j < m; j++) {
        print(entries[i][j] + " ");
      }
      println();
    }
  }
}

public Matrix IdentityMatrix(int dim) {
    Matrix id = new Matrix(dim, dim);
    for(int i = 0; i < dim; i++) {
      id.entries[i][i] = 1;
    }
    return id;
}


//*************** 3x3 linear algebra ****************
public class Matrix3x3 {
  public double a00, a01, a02;
  public double a10, a11, a12;
  public double a20, a21, a22;
  
  public Matrix3x3() {
  }

  public Matrix3x3(double... values) {
    a00 = values[0];
    a01 = values[1];
    a02 = values[2];
    a10 = values[3];
    a11 = values[4];
    a12 = values[5];
    a20 = values[6];
    a21 = values[7];
    a22 = values[8];
  }

  public double Trace() {
    return a00 + a11 + a22;
  }
  
  public Matrix3x3 Transpose() {
    return new Matrix3x3(
      a00, a10, a20,
      a01, a11, a21,
      a02, a12, a22
    );
  }

  public Matrix3x3 Copy() {
    return new Matrix3x3(
      a00,a01,a02,
      a10,a11,a12,
      a20,a21,a22
    );
  }

  public double Determinant() {
    return a00*(a11*a22-a12*a21)-a01*(a10*a22-a12*a20)+a02*(a10*a21-a11*a20);
  }

  public Matrix3x3 Inverse() {
    double det = a00*(a11*a22-a12*a21)-a01*(a10*a22-a12*a20)+a02*(a10*a21-a11*a20);
    double invD = 1 / det;
    
    return new Matrix3x3(
      invD*(a11*a22-a12*a21), invD*(a02*a21-a01*a22), invD*(a01*a12-a02*a11),
      invD*(a12*a20-a10*a22), invD*(a00*a22-a02*a20), invD*(a02*a10-a00*a12),
      invD*(a10*a21-a11*a20), invD*(a01*a20-a00*a21), invD*(a00*a11-a01*a10));
  }

  public Vector3 Solve(Vector3 b) {
    double det = a00*(a11*a22-a12*a21)-a01*(a10*a22-a12*a20)+a02*(a10*a21-a11*a20);
    double invDet = 1 / det;
    
    double b0 = b.x * invDet;
    double b1 = b.y * invDet;
    double b2 = b.z * invDet;

    return new Vector3(
      b0 * (a11*a22-a21*a12) + b1 * (a20*a12-a10*a22) + b2 * (a10*a21-a20*a11),
      b0 * (a12*a20-a10*a22) + b1 * (a00*a22-a02*a20) + b2 * (a02*a10-a00*a12),
      b0 * (a10*a21-a11*a20) + b1 * (a01*a20-a00*a21) + b2 * (a00*a11-a01*a10)
    );
  }

  public Matrix3x3 add(Matrix3x3 m) {
    return new Matrix3x3(
      a00 + m.a00, a01 + m.a01, a02 + m.a02,
      a10 + m.a10, a11 + m.a11, a12 + m.a12,
      a20 + m.a20, a21 + m.a21, a22 + m.a22
    );
  }
  
  public Matrix3x3 sub(Matrix3x3 m) {
    return new Matrix3x3(
      a00 - m.a00, a01 - m.a01, a02 - m.a02,
      a10 - m.a10, a11 - m.a11, a12 - m.a12,
      a20 - m.a20, a21 - m.a21, a22 - m.a22
    );
  }
  
  public Matrix3x3 mult(double r) {
    return new Matrix3x3(
      a00 * r, a01 * r, a02 * r,
      a10 * r, a11 * r, a12 * r,
      a20 * r, a21 * r, a22 * r
    );
  }

  public Vector3 mult(Vector3 v) {
    return new Vector3(
      v.x * a00 + v.y * a01 + v.z * a02,
      v.x * a10 + v.y * a11 + v.z * a12,
      v.x * a20 + v.y * a21 + v.z * a22
    );
  }

  public Matrix3x3 mult(Matrix3x3 m) {
    return new Matrix3x3(
      a00*m.a00 + a01*m.a10 + a02*m.a20,   a00*m.a01 + a01*m.a11 + a02*m.a21,   a00*m.a02 + a01*m.a12 + a02*m.a22,
      a10*m.a00 + a11*m.a10 + a12*m.a20,   a10*m.a01 + a11*m.a11 + a12*m.a21,   a10*m.a02 + a11*m.a12 + a12*m.a22,
      a20*m.a00 + a21*m.a10 + a22*m.a20,   a20*m.a01 + a21*m.a11 + a22*m.a21,   a20*m.a02 + a21*m.a12 + a22*m.a22);
  }

  public Matrix3x3 lin(double l1, Matrix3x3 B, double l2) {
    return new Matrix3x3(
      a00 * l1 + B.a00 * l2,
      a01 * l1 + B.a01 * l2,
      a02 * l1 + B.a02 * l2,
      a10 * l1 + B.a10 * l2,
      a11 * l1 + B.a11 * l2,
      a12 * l1 + B.a12 * l2,
      a20 * l1 + B.a20 * l2,
      a21 * l1 + B.a21 * l2,
      a22 * l1 + B.a22 * l2
    );
  }
  
  public Vector3 row0() {
    return new Vector3(a00, a01, a02);
  }

  public Vector3 row1() {
    return new Vector3(a10, a11, a12);
  }

  public Vector3 row2() {
    return new Vector3(a20, a21, a22);
  }
  
  public Vector3 column0() {
    return new Vector3(a00, a10, a20);
  }

  public Vector3 column1() {
    return new Vector3(a01, a11, a21);
  }

  public Vector3 column2() {
    return new Vector3(a02, a12, a22);
  }

  public void printMatrix() {
    println("-------------------------");
    println(a00 + " " + a01 + " " + a02);
    println(a10 + " " + a11 + " " + a12);
    println(a20 + " " + a21 + " " + a22);
    println("-------------------------");
  }
}

public class Vector3 {
  public double x, y, z;
  
  public Vector3() {
  }
  
  public Vector3(Vector3 v) {
    x = v.x;
    y = v.y;
    z = v.z;
  }
  
  public Vector3(double x, double y, double z) {
    this.x = x;
    this.y = y;
    this.z = z;
  }

  public Vector3(Matrix m) {
    x = m.entries[0][0];
    y = m.entries[1][0];
    z = m.entries[2][0];
  }

  public Vector3 Copy() {
    return new Vector3(x,y,z);
  }

  public Vector3 add(Vector3 v) {
    return new Vector3(x + v.x, y + v.y, z + v.z);
  }
  
  public Vector3 div(double r) {
    return new Vector3(x/r, y/r, z/r);
  }
  
  public void addTo(Vector3 v) {
    x+=v.x;
    y+=v.y;
    z+=v.z;
  }
  
  public void addTo(double x, double y, double z) {
    x+=x;
    y+=y;
    z+=z;
  }

  public Vector3 add(double x2, double y2, double z2) {
    return new Vector3(x + x2, y + y2, z + z2);
  }

  public Vector3 sub(Vector3 v) {
    return new Vector3(x - v.x, y - v.y, z - v.z);
  }

  public Vector3 mult(double r) {
    return new Vector3(x * r, y * r, z * r);
  }

  public double dot(Vector3 v) {
    return x * v.x + y * v.y + z * v.z;
  }

  public Vector3 cross(Vector3 v) {
    return new Vector3(
      y * v.z - z * v.y,
      z * v.x - x * v.z,
      x * v.y - y * v.x
    );
  }

  public Vector3 unit() {
    double mag = Math.sqrt(x*x + y*y + z*z);
    if(mag == 0) return new Vector3(1,0,0);
    double invMag = 1 / mag;
    return new Vector3(x*invMag, y*invMag, z*invMag);
  }

  public double mag() {
    return Math.sqrt(x*x + y*y + z*z);
  }

  public double sqr() {
    return x*x + y*y + z*z;
  }

  public Vector3 lin(double l1, Vector3 v2, double l2) {
    return new Vector3(
      x*l1 + v2.x*l2,
      y*l1 + v2.y*l2,
      z*l1 + v2.z*l2
    );
  }

  public Vector3 rotate(Vector3 axis, double angle) {
    return GetRotationMatrix_unitAxisAngle(axis.unit(), angle).mult(this);
  }

  public void printVector() {
    println("" + x + " " + y + " " + z);
  }
  
  public PVector toPVector() {
    return new PVector((float)x, (float)y, (float)z);
  }
}

public Matrix3x3 Identity3x3() {
    return new Matrix3x3(
      1, 0, 0,
      0, 1, 0,
      0, 0, 1
    );
}

public Vector3 vector(double radius, double angle) {
  return new Vector3(Math.cos(angle) * radius, Math.sin(angle) * radius, 0);
}

public Vector3 mult(Vector3 v, double a) {
  return new Vector3(v.x * a, v.y * a, v.z * a);
}

public Vector3 mult(double a, Vector3 v) {
  return new Vector3(v.x * a, v.y * a, v.z * a);
}

public Vector3 div(Vector3 v, double a) {
  return new Vector3(v.x / a, v.y / a, v.z / a);
}

public Vector3 add(Vector3 a, Vector3 b) {
  return new Vector3(a.x + b.x, a.y + b.y, a.z + b.z);
}

public Vector3 add(Vector3 a, Vector3 b, Vector3 c) {
  return new Vector3(
    a.x + b.x + c.x,
    a.y + b.y + c.y,
    a.z + b.z + c.z);
}

public Vector3 sub(Vector3 a, Vector3 b) {
  return new Vector3(a.x - b.x, a.y - b.y, a.z - b.z);
}

public Vector3 subComponent(Vector3 a, Vector3 b) {
  double magBsqr = b.x * b.x + b.y * b.y * b.z * b.z;
  if(magBsqr == 0) return a;
  double product = (a.x * b.x + a.y * b.y + a.z * b.z) / magBsqr;
  return new Vector3(a.x - product * b.x, a.y - product * b.y, a.z - product * b.z);
}

public double dist(Vector3 a, Vector3 b) {
  return sub(a, b).mag();
}

public double dot(Vector3 a, Vector3 b) {
  return a.x * b.x + a.y * b.y + a.z * b.z;
}

public Vector3 cross(Vector3 a, Vector3 b) {
  return new Vector3(
    a.y * b.z - a.z * b.y,
    a.z * b.x - a.x * b.z,
    a.x * b.y - a.y * b.x);
}

public void addTo(Vector3 a, Vector3 b, double cb) {
  a.x += b.x * cb;
  a.y += b.y * cb;
  a.z += b.z * cb;
}

public Vector3 lin(Vector3 a, double ca, Vector3 b, double cb) {
  return new Vector3(
    a.x * ca + b.x * cb,
    a.y * ca + b.y * cb,
    a.z * ca + b.z * cb);
}

public Vector3 lin(Vector3 a, double ca, Vector3 b, double cb, Vector3 c, double cc) {
  return new Vector3(
    a.x * ca + b.x * cb + c.x * cc,
    a.y * ca + b.y * cb + c.y * cc,
    a.z * ca + b.z * cb + c.z * cc);
}

public Vector3 lin(Vector3 a, double ca, Vector3 b, double cb, Vector3 c, double cc, Vector3 d, double cd) {
  return new Vector3(
    a.x * ca + b.x * cb + c.x * cc + d.x * cd,
    a.y * ca + b.y * cb + c.y * cc + d.y * cd,
    a.z * ca + b.z * cb + c.z * cc + d.z * cd);
}

public Vector3 lin(Vector3 a, double ca, Vector3 b, double cb, Vector3 c, double cc, Vector3 d, double cd, Vector3 e, double ce) {
  return new Vector3(
    a.x * ca + b.x * cb + c.x * cc + d.x * cd + e.x * ce,
    a.y * ca + b.y * cb + c.y * cc + d.y * cd + e.y * ce,
    a.z * ca + b.z * cb + c.z * cc + d.z * cd + e.z * ce);
}

public Vector3 lin(Vector3 a, double ca, Vector3 b, double cb, Vector3 c, double cc, Vector3 d, double cd, Vector3 e, double ce, Vector3 f, double cf) {
  return new Vector3(
    a.x * ca + b.x * cb + c.x * cc + d.x * cd + e.x * ce + f.x * cf,
    a.y * ca + b.y * cb + c.y * cc + d.y * cd + e.y * ce + f.y * cf,
    a.z * ca + b.z * cb + c.z * cc + d.z * cd + e.z * ce + f.z * cf);
}

public Vector3 RotateZ(Vector3 a, double arg) {
  double cos = Math.cos(arg);
  double sin = Math.sin(arg);
  return new Vector3(
    a.x * cos - a.y * sin,
    a.x * sin + a.y * cos,
    a.z);
}

public Vector3 Rotate(Vector3 v, Vector3 n, double arg) {
  double cos = Math.cos(arg);
  double sin = Math.sin(arg);
  n = n.unit();
  Vector3 vn = mult(n, dot(v, n));
  Vector3 nCrossV = cross(n, v);
  return lin(vn, 1, sub(v, vn), cos, nCrossV, sin);
}

Vector3 randomUnitSphere() {
  double x = randomGaussian();
  double y = randomGaussian();
  double z = randomGaussian();
  double rSqr = x * x + y * y + z * z;
  if(rSqr == 0) return new Vector3(1, 0, 0); 
  double mult = 1 / Math.sqrt(rSqr);
  return new Vector3(x * mult, y * mult, z * mult);
}

double constrain(double v, double min, double max) {
  if(v < min) return min;
  if(v > max) return max;
  return v;
}

double radians(double degrees) {
  return degrees * Math.PI / 180;
}

double degrees(double radians) {
  return radians * 180 / Math.PI;
}

double cosAngle(Vector3 a, Vector3 b) {
  return constrain(a.dot(b)/(a.mag() * b.mag()), -1, 1);
}

double angle(Vector3 a, Vector3 b) {
  return Math.acos(constrain(a.dot(b)/(a.mag() * b.mag()), -1, 1));
}

Vector3 map(Vector3 p0, Vector3 p1, double t0, double t1, double t) {
  if(t < t0) return p0;
  if(t > t1) return p1;
  if(t1 == t0) return p1;
  double m = (t-t0)/(t1-t0);
  return new Vector3(
    p0.x+(p1.x-p0.x)*m,
    p0.y+(p1.y-p0.y)*m,
    p0.z+(p1.z-p0.z)*m);
}

Vector3 mapFree(Vector3 p0, Vector3 p1, double t0, double t1, double t) {
  if(t1 == t0) return p1;
  double m = (t-t0)/(t1-t0);
  return new Vector3(
    p0.x+(p1.x-p0.x)*m,
    p0.y+(p1.y-p0.y)*m,
    p0.z+(p1.z-p0.z)*m);
}

public float mapConstrained(float t, float t0, float t1, float x0, float x1) {
  if(t > t1) return x1;
  if(t < t0) return x0;
  if(t0 == t1) return x1;
  return x0 + (x1-x0)*(t-t0)/(t1-t0);
}

public double mapConstrained(double t, double t0, double t1, double x0, double x1) {
  if(t > t1) return x1;
  if(t < t0) return x0;
  if(t0 == t1) return x1;
  return x0 + (x1-x0)*(t-t0)/(t1-t0);
}

public float mapFree(float t, float t0, float t1, float x0, float x1) {
  if(t0 == t1) return x1;
  return x0 + (x1-x0)*(t-t0)/(t1-t0);
}

public double mapFree(double t, double t0, double t1, double x0, double x1) {
  if(t0 == t1) return x1;
  return x0 + (x1-x0)*(t-t0)/(t1-t0);
}

public Vector3 RotateZ(Vector3 v, double cos, double sin) {
  return new Vector3(
    v.x*cos-v.y*sin,
    v.x*sin+v.y*cos,
    v.z
  );
}

public double round(double x) {
  return Math.round(x);
}



public Matrix3x3 GetRotationMatrix_unitAxisAngle(Vector3 n, double angle) {
    double cos = Math.cos(angle);
    double sin = Math.sin(angle);
    double nx = n.x;
    double ny = n.y;
    double nz = n.z;
    return new Matrix3x3(
        cos + nx*nx*(1-cos), nx*ny*(1-cos) - nz*sin, nx*nz*(1-cos) + ny*sin,
        ny*nx*(1-cos) + nz*sin, cos + ny*ny*(1-cos), ny*nz*(1-cos) - nx*sin,
        nz*nx*(1-cos) - ny*sin, nz*ny*(1-cos) + nx*sin, cos + nz*nz*(1-cos)
    );
}

public Matrix3x3 GetRotationMatrix_YawPitchRoll(double yaw, double pitch, double roll) {
    // TODO
    return Identity3x3();
}

public Matrix3x3 GetRotationMatrix_LookAt(Vector3 at, Vector3 from, double roll) {
    // TODO
    return Identity3x3();
}

public void CorrectOrthogonality(Matrix3x3 A) {
    // Gram-Schmidt process
    Vector3 e0 = new Vector3(A.a00, A.a10, A.a20);
    Vector3 e1 = new Vector3(A.a01, A.a11, A.a21);
    Vector3 e2 = new Vector3(A.a02, A.a12, A.a22);

    e0 = e0.unit();
    e1 = e1.sub(e0.mult(e0.dot(e1))).unit();
    e2 = e2.sub(e0.mult(e0.dot(e2)));
    e2 = e2.sub(e1.mult(e1.dot(e2))).unit();

    A.a00 = e0.x;
    A.a10 = e0.y;
    A.a20 = e0.z;

    A.a01 = e1.x;
    A.a11 = e1.y;
    A.a21 = e1.z;

    A.a02 = e2.x;
    A.a12 = e2.y;
    A.a22 = e2.z;
}

public double[] GetYawPitchRoll(Matrix3x3 basis) {
    CorrectOrthogonality(basis);
    Vector3 X = basis.column0();
    Vector3 Y = basis.column1();
    Vector3 Z = basis.column2();
    
    double roll = Math.atan2(-Z.x*X.y + Z.y*X.x, Y.x*X.y - Y.y*X.x);
    //Vector3 e1 = new Vector3(-X.y, X.x, 0);
    //Vector3 e2 = cross(X, e1);
    //double roll = Math.atan2(-dot(Y,e2),-dot(Y,e1));

    double yaw = Math.atan2(basis.a10, basis.a00);
    double pitch = Math.atan2(basis.a20, Math.sqrt(basis.a21*basis.a21 + basis.a22*basis.a22));
    return new double[] {yaw, pitch, roll};
}

public double[] GetEulerAngles_X1Z2X3(Matrix3x3 basis) {
    CorrectOrthogonality(basis);
    double alpha = Math.atan2(basis.a20, basis.a10);//
    double beta = Math.acos(basis.a00);//
    double gamma = Math.atan2(basis.a02, -basis.a01);
    return new double[] {alpha, beta, gamma};
}

public double[] GetEulerAngles_HeadingPitchRoll(Matrix3x3 basis) {
    CorrectOrthogonality(basis);
    Vector3 X = basis.column0();
    Vector3 Y = basis.column1();
    Vector3 Z = basis.column2();
    double pitch = Math.asin(X.z);
    double heading = Math.atan2(X.y,X.x);
    double roll = Math.atan2(Z.x, -Y.x);
    return new double[] {heading, pitch, roll};
}

public void setAffineTransform(PShape shape, Matrix3x3 basis, Vector3 offset) {
    shape.resetMatrix();

    // Rotation
    double[] angles = GetEulerAngles_X1Z2X3(basis);
    shape.rotate((float)angles[2], 1, 0, 0);
    shape.rotate((float)angles[1], 0, 0, 1);
    shape.rotate((float)angles[0], 1, 0, 0);

    shape.translate((float)offset.x, (float)offset.y, (float)offset.z);
}

public void setAffineTransform(PShape shape, Matrix3x3 basis, double scaleX, double scaleY, double scaleZ, Vector3 offset) {
    shape.resetMatrix();

    // Scaling
    shape.scale((float)scaleX, (float)scaleY, (float)scaleZ);

    // Rotation
    double[] angles = GetEulerAngles_X1Z2X3(basis);
    shape.rotate((float)angles[2], 1, 0, 0);
    shape.rotate((float)angles[1], 0, 0, 1);
    shape.rotate((float)angles[0], 1, 0, 0);

    // Translation
    shape.translate((float)offset.x, (float)offset.y, (float)offset.z);
}

public Matrix3x3 getBasisFromDirection(Vector3 direction) {
    Vector3 e1 = direction.unit();
    Vector3 e2 = new Vector3(1, 0, 0);
    if(Math.abs(e1.y) <= Math.abs(e1.x) && Math.abs(e1.y) <= Math.abs(e1.z)) e2 = new Vector3(0, 1, 0);
    if(Math.abs(e1.z) <= Math.abs(e1.x) && Math.abs(e1.z) <= Math.abs(e1.y)) e2 = new Vector3(0, 0, 1);
    e2 = e2.sub(e1.mult(e1.dot(e2))).unit();
    Vector3 e3 = e1.cross(e2);
    return new Matrix3x3(
        e1.x, e2.x, e3.x,
        e1.y, e2.y, e3.y,
        e1.z, e2.z, e3.z
    );
}

public Matrix3x3 GetInertiaMatrix_cylinder(double r, double l) {
  return new Matrix3x3(
    1 / 2.0 * r*r, 0, 0,
    0, 1 / 12.0 * (3*r*r + l*l), 0,
    0, 0, 1 / 12.0 * (3*r*r + l*l)
  );
}

public Matrix3x3 GetInertiaMatrix_T() {
  return new Matrix3x3(
    1, 0, 0,
    0, 2.9, 0,
    0, 0, 3
  );
}
