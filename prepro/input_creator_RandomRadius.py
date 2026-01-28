# Disques à rayons uniformes
#python input_creator.py nombre 100 0.5 1.0 uniform
# Disques à rayons gaussiens et les deux dernières valeurs sont optionnelles
#python input_creator.py domaine 10 8 0.5 1.0 gaussian 0.75 0.15
# Rayons constants
#python input_creator.py nombre 50 0.5 0.5 constant 0.5

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import sys
import math
from scipy.spatial.distance import cdist

def sample_radius(distribution='uniform', params=(0.5, 1.0)):
    """Échantillonne un rayon selon la distribution spécifiée"""
    if distribution == 'uniform':
        return np.random.uniform(params[0], params[1])
    elif distribution == 'gaussian':
        mean, std = params
        r = np.random.normal(mean, std)
        return max(r, 0.1 * mean)  # Évite les rayons négatifs ou trop petits
    elif distribution == 'constant':
        return params[0]
    else:
        raise ValueError(f"Distribution inconnue: {distribution}")

def can_place_disk(new_x, new_y, new_r, disks, width, height):
    """Vérifie si un nouveau disque peut être placé sans chevaucher et dans les limites"""
    # Vérifie les bords
    if new_x < new_r or new_x > width - new_r or new_y < new_r or new_y > height - new_r:
        return False
    
    # Vérifie les chevauchements avec les disques existants
    for (x, y, r) in disks:
        dist = math.sqrt((new_x - x)**2 + (new_y - y)**2)
        if dist < new_r + r - 1e-9:  # Petite marge numérique
            return False
    return True

def find_valid_position(new_r, disks, width, height, max_attempts=1000):
    """Cherche une position valide pour un nouveau disque au contact d'un disque existant"""
    if not disks:
        # Premier disque : place au centre si possible
        x = width / 2
        y = height / 2
        if can_place_disk(x, y, new_r, disks, width, height):
            return x, y
        # Sinon cherche aléatoirement
        for _ in range(max_attempts):
            x = np.random.uniform(new_r, width - new_r)
            y = np.random.uniform(new_r, height - new_r)
            if can_place_disk(x, y, new_r, disks, width, height):
                return x, y
        return None
    
    # Essaie de placer au contact d'un disque existant choisi aléatoirement
    for _ in range(max_attempts):
        # Choisit un disque existant au hasard
        idx = np.random.randint(len(disks))
        x0, y0, r0 = disks[idx]
        
        # Angle aléatoire
        theta = np.random.uniform(0, 2 * np.pi)
        
        # Distance entre les centres
        d = new_r + r0
        
        # Position candidate
        x_candidate = x0 + d * math.cos(theta)
        y_candidate = y0 + d * math.sin(theta)
        
        if can_place_disk(x_candidate, y_candidate, new_r, disks, width, height):
            return x_candidate, y_candidate
    
    return None

def pack_disks_variable(N=None, domain=None, radius_dist='uniform', radius_params=(0.5, 1.0),
                       square_side=None, max_attempts_per_disk=5000, verbose=False):
    """
    Génère un empilement de disques avec rayons variables.
    
    Args:
        N: Nombre de disques (mode carré ajustable)
        domain: Dimensions du domaine (width, height) (mode rectangulaire fixe)
        radius_dist: Distribution des rayons ('uniform', 'gaussian', 'constant')
        radius_params: Paramètres de la distribution
        square_side: Côté du carré (si spécifié en mode N)
        max_attempts_per_disk: Nombre max de tentatives par disque
        verbose: Affiche des informations
        
    Returns:
        disks: Liste de tuples (x, y, r)
        domain_size: (width, height)
    """
    
    disks = []  # Liste de tuples (x, y, r)
    
    if domain is not None:
        # Mode domaine fixe
        width, height = domain
        max_disks = max_attempts_per_disk  # Limite pour éviter boucle infinie
    else:
        # Mode N disques
        if square_side is not None:
            width = height = square_side
        else:
            # Estimation initiale de la taille
            avg_r = (radius_params[0] + radius_params[1]) / 2 if radius_dist == 'uniform' else radius_params[0]
            width = height = math.sqrt(N * math.pi * avg_r**2 / 0.8)  # Facteur de remplissage estimé
        max_disks = N
    
    placed = 0
    failed_attempts = 0
    max_total_failures = 100
    
    while placed < max_disks:
        # Échantillonne un nouveau rayon
        new_r = sample_radius(radius_dist, radius_params)
        
        # Cherche une position valide
        pos = find_valid_position(new_r, disks, width, height, max_attempts_per_disk)
        
        if pos is not None:
            x, y = pos
            disks.append((x, y, new_r))
            placed += 1
            failed_attempts = 0
            
            if verbose and placed % 100 == 0:
                print(f"Placed {placed} disks...")
        else:
            failed_attempts += 1
            if failed_attempts > max_total_failures:
                if verbose:
                    print(f"Arrêt après {failed_attempts} échecs consécutifs")
                break
    
    # En mode N, on ajuste éventuellement la taille
    if domain is None and square_side is None:
        # Trouve la boîte englobante minimale
        if disks:
            xs = [d[0] for d in disks]
            ys = [d[1] for d in disks]
            rs = [d[2] for d in disks]
            
            min_x = min(x - r for x, r in zip(xs, rs))
            max_x = max(x + r for x, r in zip(xs, rs))
            min_y = min(y - r for y, r in zip(ys, rs))
            max_y = max(y + r for y, r in zip(ys, rs))
            
            # Ajoute une petite marge
            margin = max(rs) * 0.1
            width = max_x - min_x + 2 * margin
            height = max_y - min_y + 2 * margin
            
            # Recentre les disques
            for i in range(len(disks)):
                x, y, r = disks[i]
                disks[i] = (x - min_x + margin, y - min_y + margin, r)
    
    if verbose:
        print(f"Total disks placed: {len(disks)}")
        if disks:
            avg_r = np.mean([d[2] for d in disks])
            print(f"Average radius: {avg_r:.3f}")
    
    return disks, (width, height)

def main():
    args = sys.argv[1:]
    
    if len(args) < 3:
        print("Usage:")
        print("  python input_creator.py nombre N r_min r_max [distribution] [dist_params...]")
        print("  python input_creator.py domaine width height r_min r_max [distribution] [dist_params...]")
        print("\nDistributions disponibles:")
        print("  uniform: r_min r_max (défaut)")
        print("  gaussian: mean std")
        print("  constant: radius")
        return
    
    mode = args[0]
    
    if mode == "nombre":
        N = int(args[1])
        r_min = float(args[2])
        r_max = float(args[3])
        
        # Paramètres de distribution
        if len(args) > 4:
            distribution = args[4]
            if distribution == "uniform":
                radius_params = (r_min, r_max)
            elif distribution == "gaussian":
                mean = float(args[5]) if len(args) > 5 else (r_min + r_max) / 2
                std = float(args[6]) if len(args) > 6 else (r_max - r_min) / 4
                radius_params = (mean, std)
            elif distribution == "constant":
                radius_val = float(args[5]) if len(args) > 5 else r_min
                radius_params = (radius_val,)
            else:
                print(f"Distribution inconnue: {distribution}, utilisation de 'uniform' par défaut")
                distribution = "uniform"
                radius_params = (r_min, r_max)
        else:
            distribution = "uniform"
            radius_params = (r_min, r_max)
        
        disks, (width, height) = pack_disks_variable(
            N=N,
            radius_dist=distribution,
            radius_params=radius_params,
            verbose=True
        )
        
        # Affichage
        fig, ax = plt.subplots(figsize=(8, 8))
        for (x, y, r) in disks:
            ax.add_patch(Circle((x, y), r, fill=False, edgecolor='blue'))
        ax.set_xlim(0, width)
        ax.set_ylim(0, height)
        ax.set_aspect('equal', adjustable='box')
        plt.title(f"{len(disks)} disks in domain {width:.2f}x{height:.2f}")
        plt.show()
        
    elif mode == "domaine":
        width = float(args[1])
        height = float(args[2])
        r_min = float(args[3])
        r_max = float(args[4])
        
        # Paramètres de distribution
        if len(args) > 5:
            distribution = args[5]
            if distribution == "uniform":
                radius_params = (r_min, r_max)
            elif distribution == "gaussian":
                mean = float(args[6]) if len(args) > 6 else (r_min + r_max) / 2
                std = float(args[7]) if len(args) > 7 else (r_max - r_min) / 4
                radius_params = (mean, std)
            elif distribution == "constant":
                radius_val = float(args[6]) if len(args) > 6 else r_min
                radius_params = (radius_val,)
            else:
                print(f"Distribution inconnue: {distribution}, utilisation de 'uniform' par défaut")
                distribution = "uniform"
                radius_params = (r_min, r_max)
        else:
            distribution = "uniform"
            radius_params = (r_min, r_max)
        
        disks, (W, H) = pack_disks_variable(
            domain=(width, height),
            radius_dist=distribution,
            radius_params=radius_params,
            verbose=True
        )
        
        # Affichage
        fig, ax = plt.subplots(figsize=(10, 6))
        for (x, y, r) in disks:
            ax.add_patch(Circle((x, y), r, fill=False, edgecolor='blue'))
        ax.set_xlim(0, width)
        ax.set_ylim(0, height)
        ax.set_aspect('equal', adjustable='box')
        plt.title(f"{len(disks)} disks in domain {width:.2f}x{height:.2f}")
        plt.show()
    
    # Sauvegarde dans un fichier
    if disks:
        mass = 1.0
        to_save = np.zeros((len(disks), 15))
        
        for i, (x, y, r) in enumerate(disks):
            to_save[i, 0] = x  # pos_x
            to_save[i, 1] = y  # pos_y
            to_save[i, 12] = r  # rayon
            to_save[i, 13] = 0.5 * mass * r**2  # inertie
            to_save[i, 14] = mass  # mass
        
        np.savetxt('sample_variable_radius.txt', to_save, fmt='%1.6e',
                  header="# pos_x, pos_y, vel_x, vel_y, acc_x, acc_y, pos_z, vel_z, acc_z, rot, vrot, arot, rayon, inertie, mass")
        print(f"Configuration sauvegardée dans 'sample_variable_radius.txt'")

if __name__ == "__main__":
    main()
