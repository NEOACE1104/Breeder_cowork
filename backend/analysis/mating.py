import numpy as np
import pandas as pd
import random

def calculate_fitness(pair, gebv_map, grm, trait_weights, alpha, beta):
    """
    Calculates the fitness of a single mating pair.
    Fitness = alpha * (weighted average of predicted progeny GEBV) - beta * (inbreeding coefficient)
    """
    p1, p2 = pair

    # Predicted progeny value
    predicted_gain = 0
    total_weight = sum(abs(w) for w in trait_weights.values())
    if total_weight == 0: total_weight = 1 # Avoid division by zero

    for trait, weight in trait_weights.items():
        gebv1 = gebv_map.get(p1, {}).get(trait, 0)
        gebv2 = gebv_map.get(p2, {}).get(trait, 0)
        mean_progeny_gebv = (gebv1 + gebv2) / 2
        predicted_gain += (mean_progeny_gebv * weight)

    # Normalize by total weight
    predicted_gain /= total_weight

    # Inbreeding/Coancestry from GRM
    diversity_score = 1 - grm.loc[p1, p2] # We want to maximize diversity (minimize coancestry)

    # Fitness score combines gain and diversity
    # We want to maximize gain and maximize diversity score.
    # The original prompt suggested `fitness = a*gain - b*inbreeding`.
    # Let's rephrase as `fitness = a*gain + b*diversity`
    fitness = alpha * predicted_gain + beta * diversity_score

    return fitness, predicted_gain, diversity_score

def create_initial_population(candidates, pop_size):
    """Creates an initial population of random mating pairs."""
    population = []
    for _ in range(pop_size):
        pair = tuple(random.sample(candidates, 2))
        population.append(pair)
    return population

def selection(population, fitness_scores, num_parents):
    """Selects the best individuals for the next generation."""
    fitness_scores = np.array(fitness_scores)
    # Use ranking selection to avoid premature convergence
    sorted_indices = np.argsort(fitness_scores)[::-1]

    parents = []
    for i in range(num_parents):
        parents.append(population[sorted_indices[i]])

    return parents

def crossover(parents, offspring_size, candidates):
    """Performs crossover to create new offspring."""
    offspring = []
    for _ in range(offspring_size):
        parent1 = random.choice(parents)
        parent2 = random.choice(parents)

        # Simple crossover: take one individual from each parent
        child = (parent1[0], parent2[1])
        if child[0] == child[1]: # Avoid selfing
            child = (parent1[1], parent2[0])
            if child[0] == child[1]: # Still selfing, pick a random one
                 child = (parent1[0], random.choice(list(set(candidates) - {parent1[0]})))
        offspring.append(child)
    return offspring

def mutation(offspring, mutation_rate, candidates):
    """Performs mutation on the offspring."""
    for i in range(len(offspring)):
        if random.random() < mutation_rate:
            pair = list(offspring[i])
            mutated_gene_idx = random.randint(0, 1)
            pair[mutated_gene_idx] = random.choice(list(set(candidates) - {pair[1-mutated_gene_idx]}))
            offspring[i] = tuple(pair)
    return offspring

def run_mating_ga(
    grm: pd.DataFrame,
    gebvs: pd.DataFrame, # Multi-trait GEBVs: index=sample_id, columns=traits
    request: dict
):
    """
    Runs the Genetic Algorithm for mating optimization.
    """
    # GA Parameters
    pop_size = 50
    n_generations = 100
    mutation_rate = 0.05
    num_parents_mating = 20 # Number of pairs to select for breeding

    # User inputs
    trait_weights = request['trait_weights']
    alpha = request['alpha']
    beta = request['beta']
    top_k = request['top_k']

    candidates = list(grm.index)
    if len(candidates) < 2:
        return [], {}

    # Convert GEBVs to a dictionary for faster lookup
    gebv_map = gebvs.to_dict(orient='index')

    # --- GA Loop ---
    # 1. Initialize population
    population = create_initial_population(candidates, pop_size)
    best_solution = None
    best_fitness = -np.inf

    for generation in range(n_generations):
        # 2. Calculate fitness
        fitness_data = [calculate_fitness(p, gebv_map, grm, trait_weights, alpha, beta) for p in population]
        fitness_scores = [f[0] for f in fitness_data]

        # Keep track of the best solution found so far
        max_fitness_idx = np.argmax(fitness_scores)
        if fitness_scores[max_fitness_idx] > best_fitness:
            best_fitness = fitness_scores[max_fitness_idx]
            best_solution = population[max_fitness_idx]

        # 3. Selection
        parents = selection(population, fitness_scores, num_parents_mating)

        # 4. Crossover
        offspring_crossover = crossover(parents, pop_size - len(parents), candidates)

        # 5. Mutation
        offspring_mutation = mutation(offspring_crossover, mutation_rate, candidates)

        # 6. Create new population
        population = parents + offspring_mutation

    # --- End of GA ---
    # Final population analysis
    final_fitness_data = [calculate_fitness(p, gebv_map, grm, trait_weights, alpha, beta) for p in population]

    results = []
    for i, pair in enumerate(population):
        fitness, gain, diversity = final_fitness_data[i]
        results.append({
            "parent1": pair[0],
            "parent2": pair[1],
            "predicted_gain": gain,
            "diversity_score": diversity,
            "fitness": fitness,
        })

    results_df = pd.DataFrame(results).drop_duplicates(subset=['parent1', 'parent2'])

    # Sort by fitness and take top_k
    top_pairs = results_df.sort_values(by='fitness', ascending=False).head(top_k)

    plot_data = {
        "gains": top_pairs['predicted_gain'].tolist(),
        "diversities": top_pairs['diversity_score'].tolist()
    }

    return top_pairs.to_dict('records'), plot_data
