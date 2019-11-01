using Assets.Scripts;
using System;
using System.Collections.Generic;
using System.Text;

public class DNA<T>
{
	public List<IChromosome<T>> Genes { get; private set; }
	public double Fitness { get; private set; }

	private System.Random random;
	private Func<T> getRandomGene;
	private Func<int, double> fitnessFunction;
    private Func<T, T, int> evaluationFunction;
    private Func<DNA<T>, DNA<T>,DNA<T>> crossOverBFunction;
    private Func<List<IChromosome<T>>, double> IndividualFitnessFunction;
    private readonly bool useCrossOverB;
    private float IndividualMutationRate;
    private float ChromosomeMutationRate;

	public DNA(int size, System.Random random, Func<int, double> fitnessFunction, float individualMutationRate, float chromosomeMutationRate, List<IChromosome<T>> Genes = null, bool shouldInitGenes = true, 
                bool useCrossOverB = false, Func<DNA<T>, DNA<T>, DNA<T>> crossOverBFunction = null, Func<List<IChromosome<T>>, double> individualFitnessFunction = null)
	{
        if (Genes != null)
        {
            Genes.Shuffle();
            this.Genes = new List<IChromosome<T>>(Genes);
        }
        else
            this.Genes = new List<IChromosome<T>>(size);

		this.random = random;
		this.fitnessFunction = fitnessFunction;
        IndividualMutationRate = individualMutationRate;
        ChromosomeMutationRate = chromosomeMutationRate;
        IndividualFitnessFunction = individualFitnessFunction;
        Fitness = double.NegativeInfinity;

        if (useCrossOverB)
            this.crossOverBFunction = crossOverBFunction;
	}

	public double CalculateFitness(int index)
	{
		Fitness = fitnessFunction(index);
		return Fitness;
	}

    public DNA<T> Crossover(DNA<T> otherParent)
    {
        if (useCrossOverB)
            return crossOverBFunction(this, otherParent);

        return SpliceCrossover(otherParent);
    }

	public DNA<T> SpliceCrossover(DNA<T> otherParent)
	{
		DNA<T> child = new DNA<T>(Genes.Count, random, fitnessFunction, IndividualMutationRate, ChromosomeMutationRate, 
                                    individualFitnessFunction: IndividualFitnessFunction);

        int crossOverPoint = random.Next(0, Genes.Count);
        int crossOverPoint2 = random.Next(0, Genes.Count);

        int end = Math.Max(crossOverPoint, crossOverPoint2);
        int start = Math.Min(crossOverPoint, crossOverPoint2);

        child.Genes.AddRange(Genes.GetRange(start, end - start));

		for (int i = 0; i < otherParent.Genes.Count; i++)
		{
            if (child.Genes.IndexOf(otherParent.Genes[i]) == -1)
            {
                child.Genes.Add(otherParent.Genes[i]);
            }
		}

        child.Mutate();

        return child;
	}

    public void Mutate()
    {
        if (IndividualFitnessFunction != null)
        {
            TwoOptCross();
            CalculateFiness();
            return;
        }

        for (int i = 0; i < Genes.Count; i++)
        {
            if (random.NextDouble() < ChromosomeMutationRate)
            {
                Swap(i, random.Next(0, Genes.Count));

                return;
            }
        }

        CalculateFiness();
    }

    private void TwoOptCross()
    {
        int N = Genes.Count;
        var currentRouteFitness = IndividualFitnessFunction(Genes);

        for (int i = 1; i < N - 2; i++)
        {
            for (int j = i + 1; j < N + 1; j++)
            {
                if (j - i == 1) continue;
                var newRoute = Genes.GetRange(0, i);
                var toBeReversed = Genes.GetRange(i, j - i);
                toBeReversed.Reverse();
                newRoute.AddRange(toBeReversed);
                newRoute.AddRange(Genes.GetRange(j, N - j));

                var newRouteFitness = IndividualFitnessFunction(newRoute);
                if (newRouteFitness < currentRouteFitness)
                {
                    Genes = newRoute;
                    Fitness = newRouteFitness;
                    return;
                }
            }
        }
    }

    public void AssignNewGenes(List<IChromosome<T>> Genes)
    {
        this.Genes = Genes;
    }

    public void CalculateFiness()
    {
        Fitness = IndividualFitnessFunction(Genes);
    }

    public void Swap(int indexA, int indexB)
    {
        IChromosome<T> temp = Genes[indexA];
        Genes[indexA] = Genes[indexB];
        Genes[indexB] = temp;
    }

    public override string ToString()
    {
        StringBuilder sb = new StringBuilder();

        foreach (IChromosome<T> item in Genes)
            sb.Append(item.ToString() + ", ");

        return sb.ToString();
    }
}