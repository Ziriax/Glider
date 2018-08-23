namespace FlightSim.NET
{
    public static class MathUtil
    {
        //------------------------------------------------------------------------//
        // Misc. Constants
        //------------------------------------------------------------------------//
        public const float pi = 3.14159265f;
        public const float g = -32.174f;       // acceleration due to gravity, ft/s^2
        public const float rho = 0.0023769f;   // desity of air at sea level, slugs/ft^3
        public const float tol = 0.0001f;      // float type tolerance 

        //------------------------------------------------------------------------//
        // Misc. Functions
        //------------------------------------------------------------------------//
        public static float DegreesToRadians(float deg)
        {
            return deg * pi / 180.0f;
        }

        public static float RadiansToDegrees(float rad)
        {
            return rad * 180.0f / pi;
        }
    };
}