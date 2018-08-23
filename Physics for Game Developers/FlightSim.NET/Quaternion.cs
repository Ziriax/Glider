using System;

namespace FlightSim.NET
{
    public struct Quaternion
    {
        public float n;    // number (scalar) part
        public Vector v;   // vector part: v.x, v.y, v.z

        public Quaternion(float e0, float e1, float e2, float e3)
        {
            n = e0;
            v.x = e1;
            v.y = e2;
            v.z = e3;
        }

        public float Magnitude()
        {
            return (float)Math.Sqrt(n * n + v.x * v.x + v.y * v.y + v.z * v.z);
        }

        public Vector GetVector()
        {
            return new Vector(v.x, v.y, v.z);
        }

        public float GetScalar()
        {
            return n;
        }

        // Conjugate
        public static Quaternion operator ~(Quaternion q)
        {
            return new Quaternion(q.n, -q.v.x, -q.v.y, -q.v.z);
        }

        public static Quaternion operator +(Quaternion q1, Quaternion q2)
        {
            return new Quaternion(q1.n + q2.n,
                q1.v.x + q2.v.x,
                q1.v.y + q2.v.y,
                q1.v.z + q2.v.z);
        }

        public static Quaternion operator -(Quaternion q1, Quaternion q2)
        {
            return new Quaternion(q1.n - q2.n,
                q1.v.x - q2.v.x,
                q1.v.y - q2.v.y,
                q1.v.z - q2.v.z);
        }

        public static Quaternion operator *(Quaternion q1, Quaternion q2)
        {
            return new Quaternion(q1.n * q2.n - q1.v.x * q2.v.x - q1.v.y * q2.v.y - q1.v.z * q2.v.z,
                q1.n * q2.v.x + q1.v.x * q2.n + q1.v.y * q2.v.z - q1.v.z * q2.v.y,
                q1.n * q2.v.y + q1.v.y * q2.n + q1.v.z * q2.v.x - q1.v.x * q2.v.z,
                q1.n * q2.v.z + q1.v.z * q2.n + q1.v.x * q2.v.y - q1.v.y * q2.v.x);
        }

        public static Quaternion operator *(Quaternion q, float s)
        {
            return new Quaternion(q.n * s, q.v.x * s, q.v.y * s, q.v.z * s);
        }

        public static Quaternion operator *(float s, Quaternion q)
        {
            return new Quaternion(q.n * s, q.v.x * s, q.v.y * s, q.v.z * s);
        }

        public static Quaternion operator *(Quaternion q, Vector v)
        {
            return new Quaternion(-(q.v.x * v.x + q.v.y * v.y + q.v.z * v.z),
                q.n * v.x + q.v.y * v.z - q.v.z * v.y,
                q.n * v.y + q.v.z * v.x - q.v.x * v.z,
                q.n * v.z + q.v.x * v.y - q.v.y * v.x);
        }

        public static Quaternion operator *(Vector v, Quaternion q)
        {
            return new Quaternion(-(q.v.x * v.x + q.v.y * v.y + q.v.z * v.z),
                q.n * v.x + q.v.z * v.y - q.v.y * v.z,
                q.n * v.y + q.v.x * v.z - q.v.z * v.x,
                q.n * v.z + q.v.y * v.x - q.v.x * v.y);
        }

        public static Quaternion operator /(Quaternion q, float s)
        {
            return new Quaternion(q.n / s, q.v.x / s, q.v.y / s, q.v.z / s);
        }

        static float QGetAngle(Quaternion q)
        {
            return (float)(2 * Math.Acos(q.n));
        }

        static Vector QGetAxis(Quaternion q)
        {
            Vector v;
            float m;

            v = q.GetVector();
            m = v.Magnitude();

            if(m <= MathUtil.tol)
                return new Vector();
            else
                return v / m;
        }

        static Quaternion QRotate(Quaternion q1, Quaternion q2)
        {
            return q1 * q2 * (~q1);
        }

        static Vector QVRotate(Quaternion q, Vector v)
        {
            Quaternion t;


            t = q * v * (~q);

            return t.GetVector();
        }

        static Quaternion MakeQFromEulerAngles(float x, float y, float z)
        {
            Quaternion q;
            double roll = MathUtil.DegreesToRadians(x);
            double pitch = MathUtil.DegreesToRadians(y);
            double yaw = MathUtil.DegreesToRadians(z);

            double cyaw, cpitch, croll, syaw, spitch, sroll;
            double cyawcpitch, syawspitch, cyawspitch, syawcpitch;

            cyaw = Math.Cos(0.5f * yaw);
            cpitch = Math.Cos(0.5f * pitch);
            croll = Math.Cos(0.5f * roll);
            syaw = Math.Sin(0.5f * yaw);
            spitch = Math.Sin(0.5f * pitch);
            sroll = Math.Sin(0.5f * roll);

            cyawcpitch = cyaw * cpitch;
            syawspitch = syaw * spitch;
            cyawspitch = cyaw * spitch;
            syawcpitch = syaw * cpitch;

            q.n = (float)(cyawcpitch * croll + syawspitch * sroll);
            q.v.x = (float)(cyawcpitch * sroll - syawspitch * croll);
            q.v.y = (float)(cyawspitch * croll + syawcpitch * sroll);
            q.v.z = (float)(syawcpitch * croll - cyawspitch * sroll);

            return q;
        }

        static Vector MakeEulerAnglesFromQ(Quaternion q)
        {
            double r11, r21, r31, r32, r33, r12, r13;
            double q00, q11, q22, q33;
            double tmp;
            Vector u;

            q00 = q.n * q.n;
            q11 = q.v.x * q.v.x;
            q22 = q.v.y * q.v.y;
            q33 = q.v.z * q.v.z;

            r11 = q00 + q11 - q22 - q33;
            r21 = 2 * (q.v.x * q.v.y + q.n * q.v.z);
            r31 = 2 * (q.v.x * q.v.z - q.n * q.v.y);
            r32 = 2 * (q.v.y * q.v.z + q.n * q.v.x);
            r33 = q00 - q11 - q22 + q33;

            tmp = Math.Abs(r31);
            if(tmp > 0.999999)
            {
                r12 = 2 * (q.v.x * q.v.y - q.n * q.v.z);
                r13 = 2 * (q.v.x * q.v.z + q.n * q.v.y);

                u.x = MathUtil.RadiansToDegrees(0.0f); //roll
                u.y = MathUtil.RadiansToDegrees((float)(-(MathUtil.pi / 2) * r31 / tmp)); // pitch
                u.z = MathUtil.RadiansToDegrees((float)Math.Atan2(-r12, -r31 * r13)); // yaw
                return u;
            }

            u.x = MathUtil.RadiansToDegrees((float)Math.Atan2(r32, r33)); // roll
            u.y = MathUtil.RadiansToDegrees((float)Math.Asin(-r31));         // pitch
            u.z = MathUtil.RadiansToDegrees((float)Math.Atan2(r21, r11)); // yaw
            return u;
        }

        public override string ToString()
        {
            return $"[{v.x},{v.y},{v.z},{n}]";
        }
    }
}