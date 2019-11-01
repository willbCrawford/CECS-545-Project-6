using System;

namespace Assets.Scripts
{
    public abstract class IChromosome<T>
    {
        public abstract double EvaluationFunction(IChromosome<T> Item);

        public abstract T GetT();

        public abstract int GetNode();
    }
}
