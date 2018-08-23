namespace FlightSim.NET
{
    public struct Vector
    {
        public float x;
        public float y;
        public float z;

        public Vector(float xi, float yi, float zi) : this()
        {
            x = xi;
            y = yi;
            z = zi;
        }

        public float Magnitude()
        {
            return (float) System.Math.Sqrt(x * x + y * y + z * z);
        }

        public void Normalize()
        {
            float m = (float) System.Math.Sqrt(x * x + y * y + z * z);
            if (m <= MathUtil.tol)
                m = 1;

            x /= m;
            y /= m;
            z /= m;

            if (System.Math.Abs(x) < MathUtil.tol) x = 0.0f;
            if (System.Math.Abs(y) < MathUtil.tol) y = 0.0f;
            if (System.Math.Abs(z) < MathUtil.tol) z = 0.0f;
        }

        public void Reverse()
        {
            x = -x;
            y = -y;
            z = -z;
        }

        public static Vector operator -(Vector v)
        {
            return new Vector(-v.x, -v.y, -v.z);
        }

        public static Vector operator +(Vector u, Vector v)
        {
            return new Vector(u.x + v.x, u.y + v.y, u.z + v.z);
        }

        public static Vector operator -(Vector u, Vector v)
        {
            return new Vector(u.x - v.x, u.y - v.y, u.z - v.z);
        }

        // Vector cross product (u cross v)
        public static Vector operator ^(Vector u, Vector v)
        {
            return new Vector(u.y * v.z - u.z * v.y,
                -u.x * v.z + u.z * v.x,
                u.x * v.y - u.y * v.x);
        }

        // Vector dot product
        public static float operator *(Vector u, Vector v)
        {
            return (u.x * v.x + u.y * v.y + u.z * v.z);
        }

        public static Vector operator *(float s, Vector u)
        {
            return new Vector(u.x * s, u.y * s, u.z * s);
        }

        public static Vector operator *(Vector u, float s)
        {
            return new Vector(u.x * s, u.y * s, u.z * s);
        }

        public static Vector operator /(Vector u, float s)
        {
            return new Vector(u.x / s, u.y / s, u.z / s);
        }

        // triple scalar product (u dot (v cross w))
        public static float TripleScalarProduct(Vector u, Vector v, Vector w)
        {
            return ((u.x * (v.y * w.z - v.z * w.y)) +
                    (u.y * (-v.x * w.z + v.z * w.x)) +
                    (u.z * (v.x * w.y - v.y * w.x)));
        }

        public override string ToString()
        {
            return $"({x:0.00}; {y:0.00}; {z:0.00}";
        }
    }
}