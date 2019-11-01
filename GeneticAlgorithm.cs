using Assets.Scripts;
using System;
using System.Collections.Generic;

namespace Project_4.Scripts
{
    public class GeneticAlgorithm<T>
    {
        public List<DNA<T>> Population { get; private set; }
        public int Generation { get; private set; }
        public double BestFitness { get; private set; }
        public List<IChromosome<T>> BestGenes { get; private set; }

        public int Elitism;
        public float ChromosomeMutationRate;
        public float IndividualMutationRate;

        public List<IChromosome<T>> OriginalDNA { get; private set; }
        public List<IChromosome<T>> CityDNA { get; private set; }
        private bool UseElitism;
        private List<DNA<T>> newPopulation;
        private System.Random random;
        public double fitnessSum { get; private set; }
        private int dnaSize;
        private Func<T> getRandomGene;
        private Func<int, double> fitnessFunction;
        private bool UseRoundRobin;
        private bool UseCrossOverB;
        private Func<DNA<T>, DNA<T>, DNA<T>> CrossOverBFunction;
        private Func<List<IChromosome<T>>, double> IndividualFitnessFunction;
        private float PercentCrossOver;
        private float PercentMutate;
        private SelectionType SelectionType;
        private int[,] Matrix;
        private Tuple<int, int> HighestProbabilityMatrixMember;
        private readonly int AParameter = 3;
        private readonly int BParameter = 3;
        private bool useWOC;
        public bool HasStagnated;
        private int GenerationStagnationLimit;
        private double PreviousChangedFitness;
        private int PreviousFitnessGeneration;
        private double PercentToUseWOC;

        public DNA<T> WOCMember { get; private set; }

        public GeneticAlgorithm(int populationSize, int dnaSize, System.Random random, Func<int, double> fitnessFunction, SelectionType selectionType, bool useElitism = false, 
                                int elitism = 0, float chromosomeMutationRate = 0.01f, float individualMutationRate = 0.10f, List<IChromosome<T>> Individual = null, bool useRoundRobin = true, 
                                float percentCrossOver = 0.3f, float percentMutate = 0.1f, bool useCrossOverB = false,
                                Func<DNA<T>, DNA<T>, DNA<T>> crossOverBFunction = null, Func<List<IChromosome<T>>, double> individualFitnessFunction = null, bool useWOC = false, 
                                double percentToUseWOC = 0, int GenerationStagnationLimit = 0)
        {
            Generation = 1;
            if (useElitism)
                Elitism = elitism;
            else
                Elitism = 0;
            UseElitism = useElitism;
            ChromosomeMutationRate = chromosomeMutationRate;
            IndividualMutationRate = individualMutationRate;
            Population = new List<DNA<T>>(populationSize);
            SelectionType = selectionType;

            newPopulation = new List<DNA<T>>(populationSize);
            this.random = random;
            this.dnaSize = dnaSize;
            this.fitnessFunction = fitnessFunction;
            UseRoundRobin = useRoundRobin;
            UseCrossOverB = useCrossOverB;
            CrossOverBFunction = crossOverBFunction;
            PercentCrossOver = percentCrossOver;
            PercentMutate = percentMutate;

            BestGenes = new List<IChromosome<T>>(dnaSize);
            BestFitness = double.PositiveInfinity;
            IndividualFitnessFunction = individualFitnessFunction;
            this.useWOC = useWOC;
            HasStagnated = false;
            this.GenerationStagnationLimit = GenerationStagnationLimit;
            PreviousFitnessGeneration = Generation;
            PercentToUseWOC = percentToUseWOC;

            if (Individual != null)
            {
                OriginalDNA = Individual;
                CityDNA = Individual;

                for (int i = 0; i < populationSize; i++)
                {
                    Population.Add(new DNA<T>(dnaSize, random, this.fitnessFunction, IndividualMutationRate, ChromosomeMutationRate, Genes: CityDNA, 
                                                useCrossOverB: useCrossOverB, crossOverBFunction: crossOverBFunction, 
                                                individualFitnessFunction: IndividualFitnessFunction));
                }

                if (this.useWOC)
                {
                    Matrix = new int[OriginalDNA.Count, OriginalDNA.Count];

                    for (int i = 0; i < OriginalDNA.Count; i++)
                    {
                        for (int j = 0; j < OriginalDNA.Count; j++)
                        {
                            Matrix[i, j] = 0;
                        }
                    }
                }
            }
            else
            {
                for (int i = 0; i < populationSize; i++)
                {
                    Population.Add(new DNA<T>(dnaSize, random, this.fitnessFunction, IndividualMutationRate, ChromosomeMutationRate, 
                                                useCrossOverB: useCrossOverB, crossOverBFunction: crossOverBFunction, 
                                                individualFitnessFunction: IndividualFitnessFunction));
                }
            }

            HighestProbabilityMatrixMember = new Tuple<int, int>(0, 1);
            
            if (useWOC)
                WOCMember = new DNA<T>(dnaSize, random, this.fitnessFunction, IndividualMutationRate, ChromosomeMutationRate,
                                                useCrossOverB: useCrossOverB, crossOverBFunction: crossOverBFunction,
                                                individualFitnessFunction: IndividualFitnessFunction);
        }

        private void AddToProbabilityMatrix(int populationMember)
        {
            for (int i = 0; i < Population[populationMember].Genes.Count; i++)
            {
                if (i + 1 == Population[populationMember].Genes.Count)
                {
                    break;
                }

                int j = Population[populationMember].Genes[i].GetNode()- 1;
                int k = Population[populationMember].Genes[i + 1].GetNode() - 1;

                Matrix[j, k]++;

                if (Matrix[j, k] > Matrix[HighestProbabilityMatrixMember.Item1, HighestProbabilityMatrixMember.Item2])
                {
                    HighestProbabilityMatrixMember = new Tuple<int, int>(j, k);
                }
            }
        }

        private void DetermineFinalSolutionFromWisdomOfCrowd()
        {
            List<IChromosome<T>> newRoute = new List<IChromosome<T>>();

            newRoute.AddRange(new List<IChromosome<T>>{ OriginalDNA[HighestProbabilityMatrixMember.Item1], OriginalDNA[HighestProbabilityMatrixMember.Item2]});

            int i = HighestProbabilityMatrixMember.Item2;
            while (newRoute.Count < OriginalDNA.Count)
            {
                Tuple<int, int> locallyHighestMember = new Tuple<int, int>(i, 0);
                int locallySmallestCost = Matrix[i, 0];
                for (int j = 0; j < OriginalDNA.Count; j++)
                {
                    if (Matrix[i, j] > locallySmallestCost)
                    {
                        locallyHighestMember = new Tuple<int, int>(i, j);
                        locallySmallestCost = Matrix[locallyHighestMember.Item1, locallyHighestMember.Item2];
                    }
                }

                if (newRoute.Contains(OriginalDNA[HighestProbabilityMatrixMember.Item2]))
                {
                    double locallySmallestCostGreedy = double.PositiveInfinity;
                    for (int j = 0; j < OriginalDNA.Count; j++)
                    {
                        if (newRoute.Contains(OriginalDNA[j])) continue;

                        double newSmallestCost = OriginalDNA[i].EvaluationFunction(OriginalDNA[j]);

                        if (newSmallestCost < locallySmallestCostGreedy)
                        {
                            locallyHighestMember = new Tuple<int, int>(i, j);
                            locallySmallestCostGreedy = newSmallestCost;
                        }
                    }
                }

                newRoute.AddRange(new List<IChromosome<T>> { OriginalDNA[locallyHighestMember.Item2] });

                i = locallyHighestMember.Item2;
            }

            List<IChromosome<T>> oldGenes = WOCMember.Genes;
            double oldFitness = WOCMember.Fitness;

            WOCMember.AssignNewGenes(newRoute);
            WOCMember.CalculateFiness();

            if (!double.IsNegativeInfinity(oldFitness) && oldFitness < WOCMember.Fitness)
            {
                WOCMember.AssignNewGenes(oldGenes);
                WOCMember.CalculateFiness();
            }

            if (Population[0].Fitness < WOCMember.Fitness)
            {
                WOCMember.AssignNewGenes(Population[0].Genes);
                WOCMember.CalculateFiness();
            }
        }

        public void InitializeFitness()
        {
            fitnessSum = 0;
            DNA<T> best = Population[0];

            for (int i = 0; i < Population.Count; i++)
            {
                fitnessSum += Population[i].CalculateFitness(i);

                if (Population[i].Fitness < best.Fitness)
                {
                    best = Population[i];
                }
            }

            BestFitness = best.Fitness;
            BestGenes = best.Genes;
        }

        private void PickParentsTopN(int index)
        {
            DNA<T> parent1 = Population[index];
            index++;
            DNA<T> parent2 = Population[index];

            DNA<T> child = parent1.Crossover(parent2);

            newPopulation.Add(child);
        }

        private List<DNA<T>> ChooseMatingPool()
        {
            List<DNA<T>> matingPool = new List<DNA<T>>(Population.Count);

            for (int i = 0; i < Population.Count; i++)
            {
                 int fitnessProportion = (int)(Population.Count * (1 - (Population[i].Fitness / fitnessSum)));
                
                for (int j = 0; j < fitnessProportion; j++)
                {
                    matingPool.Add(Population[i]);
                }
            }

            return matingPool;
        }
        
        private void BeginMating(List<DNA<T>> matingPool)
        {
            int elementA = random.Next(matingPool.Count);
            int elementB = random.Next(matingPool.Count);

            DNA<T> parentA = matingPool[elementA];
            DNA<T> parentB = matingPool[elementB];

            DNA<T> child = parentA.Crossover(parentB);

            newPopulation.Add(child);
        }

        public void ResetMatrix()
        {
            for (int i = 0; i < OriginalDNA.Count; i++)
            {
                for (int j = 0; j < OriginalDNA.Count; j++)
                {
                    Matrix[i, j] = 0;
                }
            }
        }

        public void NewGeneration(int numNewDNA = 0)
        {
            int finalCount = Population.Count + numNewDNA;

            if (finalCount <= 0)
            {
                return;
            }

            newPopulation.Clear();

            List<DNA<T>> matingPool = new List<DNA<T>>(Population.Count);
            if (SelectionType == SelectionType.ROULETTE)
            {
                matingPool = ChooseMatingPool();
                matingPool.Shuffle();
            }
            else
            {
                matingPool = null;
            }

            for (int i = 0; i < Population.Count; i++)
            {
                if (i < Elitism)
                {
                    Population[i].Mutate();
                    newPopulation.Add(Population[i]);
                }

                if (i < (int)(Population.Count * PercentCrossOver))
                {
                    if (SelectionType == SelectionType.TOP)
                        PickParentsTopN(i);
                }
                else
                {
                    if (SelectionType == SelectionType.TOP)
                    {
                        Population[i].Mutate();
                        newPopulation.Add(Population[i]);
                    }
                }

                if (SelectionType == SelectionType.ROULETTE)
                    BeginMating(matingPool);

                if (i < (int)(Population.Count * PercentToUseWOC))
                    if (useWOC)
                        AddToProbabilityMatrix(i);

                if (newPopulation.Count >= Population.Count)
                {
                    break;
                }
            }

            if (useWOC)
            {
                DetermineFinalSolutionFromWisdomOfCrowd();
                ResetMatrix();
                Population[Population.Count - 1] = WOCMember;
            }

            List<DNA<T>> tmpList = Population;
            Population = newPopulation;
            newPopulation = tmpList;

            CalculateFitness();
            Population.Sort(CompareDNA);

            Generation++;

            if (PreviousChangedFitness == BestFitness)
            {
                if (GenerationStagnationLimit == Generation - PreviousFitnessGeneration)
                {
                    HasStagnated = true;
                }
            }
            else
            {
                PreviousChangedFitness = BestFitness;
                PreviousFitnessGeneration = Generation;
            }
        }

        private int CompareDNA(DNA<T> a, DNA<T> b)
        {
            if (a.Fitness > b.Fitness)
            {
                return 1;
            }
            else if (a.Fitness < b.Fitness)
            {
                return -1;
            }
            else
            {
                return 0;
            }
        }

        private void CalculateFitness()
        {
            fitnessSum = 0;
            DNA<T> best = Population[0];

            for (int i = 0; i < Population.Count; i++)
            {
                fitnessSum += Population[i].CalculateFitness(i);

                if (Population[i].Fitness < best.Fitness)
                {
                    best = Population[i];
                }
            }

            BestFitness = best.Fitness;
            BestGenes = best.Genes;
        }

        private DNA<T> ChooseParent()
        {
            return ChooseParentRandom();
        }

        private DNA<T> ChooseParentRandom()
        {
            //double randomNumber = random.NextDouble() * fitnessSum;
            double selectedValue = random.NextDouble() * fitnessSum;
            double reciprocal = 1 / selectedValue;
            double currentSum = 0;

            for (int i = 0; i < Population.Count; i++)
            {
                currentSum += (1 / Population[i].Fitness);
                if (currentSum >= reciprocal)
                {
                    return Population[i];
                }
            }

            return Population[random.Next(Population.Count)];
        }

        private DNA<T> ChooseParentRoundRobin()
        {
            DNA<T> parent1 = Population[0];
            DNA<T> parent2 = Population[1];

            foreach (DNA<T> Individual in Population)
            {
                if (parent1.Fitness < Individual.Fitness)
                    parent1 = Individual;
                else if (parent2.Fitness < Individual.Fitness)
                    parent2 = Individual;
            }

            return parent1.Crossover(parent2);
        }
    }

    public enum SelectionType
    {
        ROULETTE,
        TOURNAMENT,
        TOP
    }
}