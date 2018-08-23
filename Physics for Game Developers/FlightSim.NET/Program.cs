using System;

namespace FlightSim.NET
{
    class Program
    {
        static void Main(string[] args)
        {
            Console.WriteLine("C#");

            const float dt = 0.1f;

            Glider glider = new Glider();

            for (int i = 0; i < 10; ++i)
            {
                Console.WriteLine($"pos: {glider.Airplane.vPosition} rot: {glider.Airplane.qOrientation}");
                glider.StepSimulation(dt);

                switch (i)
                {
                    case 5:
                        glider.PitchDown();
                        break;
                }
            }

            Console.ReadLine();
        }
    }
}
