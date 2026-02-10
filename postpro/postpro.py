import re
import numpy as np
from typing import Dict, Any, List
import os
from math import pi
from collections import defaultdict
import heapq  # Pour l'algorithme de Dijkstra

FLOAT_RE = re.compile(r'[-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?')

def read_mfloe_conf(path: str) -> Dict[str, Any]:
    floe_cols = ['pos_x','pos_y','zpos','vel_x','vel_y','zvel','acc_x','acc_y','zacc',
                 'rot','vrot','arot','radius','height','inertia','mass']
    inter_cols = ['i','j','isBonded','fn','fnb','ft','ftb','fs','fsb','A','coverage','dn0']

    meta = {}
    floe_list = []
    inter_list = []

    with open(path, 'r') as f:
        lines = [ln.rstrip('\n') for ln in f]

    i = 0
    L = len(lines)
    # helper to parse floats in a line
    def parse_floats(line):
        return [float(x) for x in FLOAT_RE.findall(line)]

    while i < L:
        line = lines[i].strip()
        if not line:
            i += 1
            continue
        parts = line.split()
        key = parts[0]

        # header MFloe version
        if key == 'MFloe':
            if len(parts) >= 2:
                meta['prog'] = 'MFloe'
                meta['version'] = parts[1]
            i += 1
            continue

        # known scalar keys (store as float when possible)
        # list des clés attendues d'après saveConf
        scalar_keys = {'t','tmax','dt','interClose','interOut','interHist','dVerlet',
                       'zgravNorm','iconf','kn','kt','mu','Gc','nDriven','FloeElements','Interactions'}
        if key in scalar_keys and len(parts) >= 2:
            # pour FloeElements et Interactions on gère ci-dessous
            if key not in ('FloeElements','Interactions'):
                # try float, otherwise store raw
                try:
                    meta[key] = float(parts[1])
                except ValueError:
                    meta[key] = parts[1]
                i += 1
                continue

        # FloeElements block
        if key == 'FloeElements':
            # format: "FloeElements N"
            try:
                n = int(parts[1])
            except Exception:
                n = 0
            i += 1
            for _ in range(n):
                if i >= L:
                    break
                ln = lines[i].strip()
                # extract all numbers from the line
                vals = parse_floats(ln)
                if len(vals) == 0:
                    # skip empty line and continue reading
                    i += 1
                    continue
                # Expect 16 values per floe (pos x,y ; zpos ; vel x,y ; zvel ; acc x,y ; zacc ;
                # rot ; vrot ; arot ; radius ; height ; inertia ; mass) => 16
                if len(vals) < 16:
                    # tolérance : laisser la ligne mais compléter par NaN
                    vals = vals + [np.nan] * (16 - len(vals))
                elif len(vals) > 16:
                    # si plus de 16 valeurs (rare), tronquer à 16
                    vals = vals[:16]
                floe_list.append(vals)
                i += 1
            continue

        # Interactions block
        if key == 'Interactions':
            try:
                m = int(parts[1])
            except Exception:
                m = 0
            i += 1
            for _ in range(m):
                if i >= L:
                    break
                ln = lines[i].strip()
                vals = parse_floats(ln)
                if len(vals) == 0:
                    i += 1
                    continue
                # Expect 13 values: i j isBonded fn fnb ft ftb fs fsb A coverage dn0
                if len(vals) < 12:
                    vals = vals + [np.nan] * (12 - len(vals))
                elif len(vals) > 12:
                    vals = vals[:12]
                inter_list.append(vals)
                i += 1
            continue

        # lignes inconnues (ex: Drivings.* écrit par driving->write). On skip jusqu'à la prochaine ligne.
        # Simplement on avance d'une ligne — ceci laisse le parseur tomber sur "FloeElements" quand il arrive.
        i += 1

    # conversion en numpy arrays
    floe_arr = np.array(floe_list, dtype=np.float64) if floe_list else np.empty((0, len(floe_cols)))
    inter_arr = np.array(inter_list, dtype=np.float64) if inter_list else np.empty((0, len(inter_cols)))

    # pour interactions, convertir colonnes i,j,isBonded en int si possible
    if inter_arr.size:
        try:
            inter_ints = inter_arr[:, :3].astype(np.int64)
            inter_arr[:, :3] = inter_ints
        except Exception:
            # si conversion échoue on laisse en float
            pass

    return {
        'meta': meta,
        'floe': floe_arr,
        'floe_cols': floe_cols,
        'interactions': inter_arr,
        'inter_cols': inter_cols
    }


def compute_systole_and_volume(pos, idx_i, idx_j, contact_dist):
    """
    Calcule la systole (longueur de la plus courte boucle non triviale) 
    et le volume (somme des longueurs des arêtes) du graphe des interactions.
    
    Args:
        pos: array (N,2) des positions des particules
        idx_i, idx_j: arrays (M,) des indices des particules en interaction
        contact_dist: array (M,) des distances entre particules en interaction
        
    Returns:
        tuple: (systole, volume) où systole est la longueur de la plus courte boucle,
               ou np.nan si aucune boucle n'existe
    """
    N = len(pos)
    M = len(idx_i)
    
    if M == 0:
        return np.nan, 0.0
    
    # Volume = somme des longueurs des arêtes
    volume = np.sum(contact_dist)
    
    # Construction du graphe comme liste d'adjacence
    graph = defaultdict(list)
    for k in range(M):
        i = int(idx_i[k])
        j = int(idx_j[k])
        d = contact_dist[k]
        graph[i].append((j, d))
        graph[j].append((i, d))
    
    # Recherche de la plus courte boucle avec BFS adapté
    min_cycle = float('inf')
    
    for start_node in graph:
        # Dictionnaire des distances depuis start_node
        dist = {start_node: 0}
        # Dictionnaire des parents pour reconstruction du chemin
        parent = {start_node: -1}
        # File pour BFS
        queue = [(0, start_node)]  # (distance, node)
        
        while queue:
            # Utilisation d'un tas pour obtenir toujours le nœud le plus proche
            current_dist, current = heapq.heappop(queue)
            
            for neighbor, weight in graph[current]:
                # Éviter de revenir directement en arrière
                if neighbor == parent[current]:
                    continue
                    
                new_dist = current_dist + weight
                
                if neighbor not in dist:
                    # Nouveau nœud découvert
                    dist[neighbor] = new_dist
                    parent[neighbor] = current
                    heapq.heappush(queue, (new_dist, neighbor))
                else:
                    # Cycle détecté !
                    cycle_length = new_dist + dist[neighbor]
                    if cycle_length < min_cycle:
                        min_cycle = cycle_length
    
    systole = min_cycle if min_cycle < float('inf') else np.nan
    
    return systole, volume

def compute_metrics(data, domain_area=None, dr=0.5, r_max=None, compute_gr=False):
    """
    Calcule un grand nombre d'estimateurs à partir de `data` (format read_mfloe_conf).
    Retour: dict contenant les métriques + tableaux utiles.

    Args:
      data: dict renvoyé par read_mfloe_conf
      domain_area: float (optionnel). Si None => use bounding box area of positions.
      dr: bin width for g(r)
      r_max: max radius for g(r). If None, set to half of min(box dim).
      compute_gr: bool, si True calcule g(r) (peut être coûteux pour grands N).
    """
    floe = data['floe']            # N x 16
    fcols = data['floe_cols']
    inter = data['interactions']   # M x 13
    icols = data['inter_cols']

    # helpers pour colonnes
    def col_idx(name):
        return fcols.index(name)
    pos = floe[:, [col_idx('pos_x'), col_idx('pos_y')]] if len(fcols)>0 else floe[:, :2]
    rads = floe[:, fcols.index('radius')]
    heights = floe[:, fcols.index('height')]
    masses = floe[:, fcols.index('mass')]
    inertias = floe[:, fcols.index('inertia')]
    vel = floe[:, [fcols.index('vel_x'), fcols.index('vel_y')]]
    vrot = floe[:, fcols.index('vrot')]

    N = pos.shape[0]
    M = inter.shape[0]

    # Convertir i,j,isBonded en int/boolean
    if M > 0:
        idx_i = inter[:, icols.index('i')].astype(int)
        idx_j = inter[:, icols.index('j')].astype(int)
        isBonded = inter[:, icols.index('isBonded')].astype(int)
        fn = inter[:, icols.index('fn')]
        ft = inter[:, icols.index('ft')]
        # autres colonnes utilisables : fnb, ftb, fs, fsb, A, coverage, Gc, dn0 ...
    else:
        idx_i = np.array([], dtype=int)
        idx_j = np.array([], dtype=int)
        isBonded = np.array([], dtype=int)
        fn = np.array([])
        ft = np.array([])

    # Domain area
    if domain_area is None:
        xmin, ymin = pos.min(axis=0)
        xmax, ymax = pos.max(axis=0)
        domain_area = max(1e-12, (xmax - xmin) * (ymax - ymin))
        box_dims = (xmax - xmin, ymax - ymin)
    else:
        box_dims = None

    # Packing fraction
    areas = pi * rads**2
    A_particles = areas.sum()
    packing_fraction = A_particles / domain_area

    # Basic stats
    stats = {}
    stats['N'] = N
    stats['M'] = M
    stats['Nbonds'] = int((isBonded == 1).sum())
    stats['fraction_bonded_of_contacts'] = (isBonded == 1).sum() / max(1, M)
    stats['fraction_bonded_of_particles'] = (isBonded == 1).sum() / max(1, N)
    stats['mean_radius'] = float(rads.mean() if N>0 else np.nan)
    stats['std_radius'] = float(rads.std() if N>0 else np.nan)
    stats['median_radius'] = float(np.median(rads) if N>0 else np.nan)
    stats['min_radius'] = float(rads.min() if N>0 else np.nan)
    stats['max_radius'] = float(rads.max() if N>0 else np.nan)
    stats['A_particles'] = float(A_particles)
    stats['domain_area'] = float(domain_area)
    stats['packing_fraction'] = float(packing_fraction)
    stats['box_dims'] = box_dims

    # Center of mass
    if N > 0:
        total_mass = masses.sum()
        com = (pos * masses[:,None]).sum(axis=0) / total_mass if total_mass>0 else pos.mean(axis=0)
    else:
        com = np.array([np.nan, np.nan])
    stats['COM'] = com

    # Kinetic energies
    if N>0:
        speeds2 = (vel**2).sum(axis=1)
        K_trans = 0.5 * (masses * speeds2).sum()
        K_rot = 0.5 * (inertias * (vrot**2)).sum()
    else:
        K_trans = K_rot = 0.0
    stats['K_trans'] = float(K_trans)
    stats['K_rot'] = float(K_rot)
    stats['K_total'] = float(K_trans + K_rot)

    # Granular temperature (variance of velocity components)
    if N>0:
        v_mean = vel.mean(axis=0)
        Tg = ((vel - v_mean)**2).sum() / (2 * N)   # moyenne par degré de liberté (2D)
    else:
        Tg = np.nan
    stats['granular_temperature'] = float(Tg)

    # Coordination number & degree distribution
    deg = np.zeros(N, dtype=int)
    for a,b in zip(idx_i, idx_j):
        if 0 <= a < N: deg[a] += 1
        if 0 <= b < N: deg[b] += 1
    if N>0:
        z_mean = deg.mean()
    else:
        z_mean = np.nan
    stats['coordination_number_mean'] = float(z_mean)
    stats['degree_hist'] = np.bincount(deg)

    # Per-contact geometry vectors and forces
    forces_on_i = np.zeros((N,2), dtype=float)   # accumulate force vector on particle i
    forces_on_j = np.zeros((N,2), dtype=float)
    contact_forces = np.zeros((M,2), dtype=float)
    contact_normal = np.zeros((M,2), dtype=float)
    contact_dist = np.zeros(M, dtype=float)
    for k in range(M):
        i = idx_i[k]; j = idx_j[k]
        if not (0 <= i < N and 0 <= j < N):
            continue
        rij = pos[j] - pos[i]          # vector i->j
        dist = np.linalg.norm(rij)
        if dist == 0:
            # fallback: tiny random vector
            n = np.array([1.0, 0.0])
        else:
            n = rij / dist
        t = np.array([-n[1], n[0]])    # in-plane perpendicular (90°)
        # assume fn is signed normal force (positive compression), ft is signed tangential
        fn_k = fn[k] if len(fn)>k else 0.0
        ft_k = ft[k] if len(ft)>k else 0.0
        f_vec = fn_k * n + ft_k * t
        contact_forces[k] = f_vec
        contact_normal[k] = n
        contact_dist[k] = dist
        # force on i from j is +f_vec, on j from i is -f_vec
        forces_on_i[i] += f_vec
        forces_on_j[j] += -f_vec

    net_forces = forces_on_i + forces_on_j*0  # net forces per particle (forces_on_j is negative on j only)
    stats['max_net_force'] = float(np.linalg.norm(net_forces, axis=1).max() if N>0 else 0.0)
    stats['mean_net_force'] = float(np.linalg.norm(net_forces, axis=1).mean() if N>0 else 0.0)

    # Calcul de la systole et du volume du graphe
    systole, graph_volume = compute_systole_and_volume(pos, idx_i, idx_j, contact_dist)
    
    # Calcul de la systolicité
    systolicity = graph_volume / systole if not np.isnan(systole) and systole > 0 else np.nan

    # Ajouter ces métriques aux stats
    stats['graph_volume'] = float(graph_volume)
    stats['systole'] = float(systole)
    stats['systolicity'] = float(systolicity)


    # Contact force statistics
    fn_all = np.array(fn) if M>0 else np.array([])
    ft_all = np.array(ft) if M>0 else np.array([])
    fmag = np.linalg.norm(contact_forces, axis=1) if M>0 else np.array([])
    stats['fn_mean'] = float(np.mean(fn_all) if fn_all.size else np.nan)
    stats['ft_mean'] = float(np.mean(ft_all) if ft_all.size else np.nan)
    stats['fmag_mean'] = float(np.mean(fmag) if fmag.size else np.nan)
    stats['fn_hist'] = None
    stats['ft_hist'] = None

    # Tenseur des contraintes (virial) 2x2 : sigma = (1/Area) sum(r_ij ⊗ f_ij)
    sigma = np.zeros((2,2))
    if M>0:
        for k in range(M):
            rvec = pos[idx_j[k]] - pos[idx_i[k]]       # r_ij
            fvec = contact_forces[k]                   # force on i from j
            sigma += np.outer(rvec, fvec)
        sigma = sigma / domain_area
    stats['stress_tensor'] = sigma
    stats['pressure_estimate'] = -0.5 * np.trace(sigma)  # sign convention: positive pressure if compressive

    # Fabric tensor (avg n ⊗ n) and anisotropy
    fabric = np.zeros((2,2))
    if M>0:
        for k in range(M):
            n = contact_normal[k]
            fabric += np.outer(n, n)
        fabric = fabric / M
        # anisotropy = difference of eigenvalues
        evals = np.linalg.eigvalsh(fabric)
        anisotropy = float(abs(evals[1] - evals[0]))
    else:
        anisotropy = np.nan
    stats['fabric'] = fabric
    stats['fabric_anisotropy'] = anisotropy

    # Cluster analysis based on bonds (isBonded==1)
    parent = np.arange(N, dtype=int)
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    def union(a,b):
        ra = find(a); rb = find(b)
        if ra != rb:
            parent[rb] = ra

    for k in range(M):
        if isBonded[k] == 1:
            a = idx_i[k]; b = idx_j[k]
            if 0 <= a < N and 0 <= b < N:
                union(a,b)
    roots = [find(i) for i in range(N)]
    from collections import Counter
    cluster_counts = Counter(roots)
    cluster_sizes = np.array(list(cluster_counts.values())) if len(cluster_counts)>0 else np.array([])
    if cluster_sizes.size:
        largest_cluster = cluster_sizes.max()
        mean_cluster_size = cluster_sizes.mean()
    else:
        largest_cluster = 0
        mean_cluster_size = 0
    stats['n_clusters'] = int(len(cluster_sizes))
    stats['largest_cluster_size'] = int(largest_cluster)
    stats['mean_cluster_size'] = float(mean_cluster_size)
    stats['cluster_size_hist'] = None

    # Radial distribution function g(r)
    gr = None
    gr_r = None
    if compute_gr and N>1:
        # compute all pairwise distances (O(N^2) !)
        # choose r_max if not provided
        if r_max is None and box_dims is not None:
            r_max = 0.5 * min(box_dims)
        elif r_max is None:
            r_max = np.linalg.norm(pos.max(axis=0)-pos.min(axis=0)) / 2.0
        nbins = max(1, int(np.ceil(r_max / dr)))
        edges = np.linspace(0, r_max, nbins+1)
        counts = np.zeros(nbins, dtype=int)
        # compute pairwise distances efficiently
        for i in range(N):
            d = np.linalg.norm(pos[i+1:] - pos[i], axis=1)
            # only keep <= r_max
            d = d[d>0]
            idxs = (d / dr).astype(int)
            idxs = idxs[idxs < nbins]
            for idxx in idxs:
                counts[idxx] += 1
        # normalization: g(r) = counts / (2π r dr * rho * N)
        rho = N / domain_area
        r_centers = 0.5 * (edges[:-1] + edges[1:])
        shell_areas = 2 * pi * r_centers * dr
        norm = rho * N * shell_areas
        with np.errstate(divide='ignore', invalid='ignore'):
            g_r = counts / norm
        gr = g_r
        gr_r = r_centers
        stats['g_r_r'] = gr_r
        stats['g_r'] = gr

    # summary histograms
    stats['radius_hist'] = np.histogram(rads, bins=20)
    stats['force_magnitude_hist'] = np.histogram(fmag, bins=30) if fmag.size else (np.array([]), np.array([]))
    stats['cluster_sizes'] = cluster_sizes

    # package results
    results = {
        'stats': stats,
        'net_forces': net_forces,
        'contact_forces': contact_forces,
        'contact_normals': contact_normal,
        'contact_distances': contact_dist,
        'degree_array': deg,
        'cluster_sizes': cluster_sizes,
        'g_r': (gr_r, gr) if gr is not None else None
    }
    return results




# --- suppose que read_mfloe_conf et compute_metrics sont déjà définies/importées ---
# from ton_module import read_mfloe_conf, compute_metrics

def save_concat_metrics(
    liste_conf: List[str],
    out_prefix: str = "summary",
    dr: float = 0.2,
    r_max: float = None,
    domain_area: float = None,
    verbose: bool = True
):
    """
    Parcourt la liste de fichiers conf, calcule metrics (avec g(r)) et sauvegarde :
      - out_prefix + '_summary.csv' : résumé (une ligne par fichier)
      - out_prefix + '_g_r_matrix.txt' : matrice g(r) (colonnes: [r_center, g_r_file1, g_r_file2, ...])
    Choix: calcule un r_max global si non fourni pour garantir des bins identiques pour tous les fichiers.
    """
    if len(liste_conf) == 0:
        raise ValueError("liste_conf vide")

    # 1) estimation du domaine global (pour r_max cohérent) : on lit juste positions pour chaque conf
    xmin = ymin = np.inf
    xmax = ymax = -np.inf
    for file in liste_conf:
        data = read_mfloe_conf(file)
        floe = data['floe']
        if floe.shape[0] == 0:
            continue
        # indices pos_x,pos_y
        pos_x = floe[:, data['floe_cols'].index('pos_x')]
        pos_y = floe[:, data['floe_cols'].index('pos_y')]
        xmin = min(xmin, pos_x.min())
        ymin = min(ymin, pos_y.min())
        xmax = max(xmax, pos_x.max())
        ymax = max(ymax, pos_y.max())

    if xmin == np.inf:
        raise RuntimeError("Aucune particule trouvée dans les fichiers fournis")

    box_x = xmax - xmin
    box_y = ymax - ymin
    if r_max is None:
        r_max = 0.5 * min(box_x, box_y)
        if verbose:
            print(f"Choix automatique r_max = {r_max:.4g} (half min box dims {box_x:.4g} x {box_y:.4g})")

    # compute number of bins and r_centers (same pour tous)
    nbins = max(1, int(np.ceil(r_max / dr)))
    edges = np.linspace(0, r_max, nbins + 1)
    r_centers = 0.5 * (edges[:-1] + edges[1:])

    # Prepare output containers
    scalar_rows = []
    gr_cols = []  # each element: g_r array (length nbins)
    filenames_for_rows = []

    # header / order for summary columns
    summary_keys = [
    'file', 't', 'N', 'M', 'Nbonds',
    'fraction_bonded_of_contacts', 'fraction_bonded_of_particles',
    'mean_radius', 'std_radius', 'packing_fraction',
    'COM_x', 'COM_y',
    'K_trans', 'K_rot', 'K_total',
    'granular_temperature', 'coordination_number_mean',
    'max_net_force', 'mean_net_force',
    'fn_mean', 'ft_mean', 'fmag_mean',
    'pressure_estimate', 'fabric_anisotropy',
    'n_clusters', 'largest_cluster_size', 'mean_cluster_size',
    'graph_volume', 'systole', 'systolicity'
    ]

    # 2) boucle principale — calcul des metrics avec compute_gr True et r_max fixé
    for file in liste_conf:
        if verbose:
            print("Processing:", file)
        data = read_mfloe_conf(file)
        # Time if available
        t_val = data.get('meta', {}).get('t', np.nan)

        # compute metrics with same r_max and dr
        res = compute_metrics(data, domain_area=domain_area, dr=dr, r_max=r_max, compute_gr=True)

        stats = res['stats']
        # extract scalars (avec fallback np.nan)
        row = [
            os.path.basename(file),                        # file
            stats.get('t', t_val) if 't' in stats else t_val,  # t (try stats then meta)
            stats.get('N', np.nan),
            stats.get('M', np.nan),
            stats.get('Nbonds', np.nan),
            stats.get('fraction_bonded_of_contacts', np.nan),
            stats.get('fraction_bonded_of_particles', np.nan),
            stats.get('mean_radius', np.nan),
            stats.get('std_radius', np.nan),
            stats.get('packing_fraction', np.nan),
            stats.get('COM')[0] if 'COM' in stats and stats['COM'] is not None else np.nan,
            stats.get('COM')[1] if 'COM' in stats and stats['COM'] is not None else np.nan,
            stats.get('K_trans', np.nan),
            stats.get('K_rot', np.nan),
            stats.get('K_total', np.nan),
            stats.get('granular_temperature', np.nan),
            stats.get('coordination_number_mean', np.nan),
            stats.get('max_net_force', np.nan),
            stats.get('mean_net_force', np.nan),
            stats.get('fn_mean', np.nan),
            stats.get('ft_mean', np.nan),
            stats.get('fmag_mean', np.nan),
            stats.get('pressure_estimate', np.nan),
            stats.get('fabric_anisotropy', np.nan),
            stats.get('n_clusters', np.nan),
            stats.get('largest_cluster_size', np.nan),
            stats.get('mean_cluster_size', np.nan),
            stats.get('graph_volume', np.nan),       # <-- Ajouter
            stats.get('systole', np.nan),           # <-- Ajouter
            stats.get('systolicity', np.nan)        # <-- Ajouter
        ]
        scalar_rows.append(row)
        filenames_for_rows.append(os.path.basename(file))

        # g(r) : res['g_r'] est (r_centers, g_r) ou None
        gr_pair = res.get('g_r', None)
        if gr_pair is None:
            # mettre NaNs
            gr_cols.append(np.full(nbins, np.nan))
        else:
            gr_r_centers, gr_vals = gr_pair
            # gr_vals peut être plus court si r_max trop grand; on aligne sur nbins :
            if len(gr_vals) == nbins:
                gr_cols.append(gr_vals)
            else:
                # si compute_metrics a utilisé exactement les mêmes r_max/dr on est OK,
                # sinon on interpole sur r_centers global
                gr_cols.append(np.interp(r_centers, gr_r_centers, gr_vals, left=np.nan, right=np.nan))

    # 3) sauvegarde du résumé (CSV) via np.savetxt ; on convertit tout en chaînes pour garder le nom de fichier
    summary_header = ",".join(summary_keys)
    # build array-of-strings
    str_rows = []
    for row in scalar_rows:
        # str conversion with reasonable precision
        str_row = [str(row[0])] + [("{:.8e}".format(x) if (isinstance(x, (int, float, np.floating, np.integer)) and not np.isnan(x)) else "nan") for x in row[1:]]
        str_rows.append(str_row)
    arr_str = np.array(str_rows, dtype=object)
    summary_path = out_prefix + "_summary.csv"
    np.savetxt(summary_path, arr_str, delimiter=",", fmt="%s", header=summary_header, comments='')
    if verbose:
        print("Saved summary to", summary_path)

    # 4) sauvegarde de g(r) matrix : colonnes = [r_centers, g_r_file1, g_r_file2, ...]
    gr_matrix = np.vstack(gr_cols).T  # shape (nbins, nfiles)
    # construct 2D array with first column = r_centers
    out_mat = np.column_stack([r_centers, gr_matrix])
    # header for g(r) matrix
    gr_header = "r_center," + ",".join([os.path.basename(f) for f in liste_conf])
    gr_path = out_prefix + "_g_r_matrix.txt"
    np.savetxt(gr_path, out_mat, delimiter=",", header=gr_header, comments='', fmt='%.6e')
    if verbose:
        print("Saved g(r) matrix to", gr_path)

    return summary_path, gr_path

# --- Exemple d'utilisation ---
# liste_conf = ["run_t000.conf", "run_t010.conf", "run_t020.conf"]
# save_concat_metrics(liste_conf, out_prefix="sim_all", dr=0.2)


def creer_dossier(nom_dossier):
    try:
        # Verifiez si le dossier n'existe pas deja
        if not os.path.exists(nom_dossier):
            # Creez le dossier
            os.makedirs(nom_dossier)
            print(f"Dossier '{nom_dossier}' cree avec succes.")
        else:
            print(f"Le dossier '{nom_dossier}' existe deja.")
    except Exception as e:
        print(f"Erreur lors de la creation du dossier '{nom_dossier}': {str(e)}")

def num_key(name):
    m = re.search(r'\d+', name)
    return int(m.group()) if m else -1   # -1 pour mettre les noms sans chiffre devant

nom_dossier_simulations = "./surfaceenergy_frottement/"
liste_simulation = os.listdir(nom_dossier_simulations)
liste_simulation = [f for f in liste_simulation if "surface" in f]


for i in liste_simulation:
    dossier_to_analyze = nom_dossier_simulations + i
    liste_conf = sorted(os.listdir(dossier_to_analyze), key = num_key)
    liste_conf = [dossier_to_analyze + "/" + f for f in liste_conf if "conf" in f]
    
    ###
    #liste_conf = liste_conf[:2]
    ###
    creer_dossier("./postprocess/" + i)
    save_concat_metrics(liste_conf, out_prefix="./postprocess/" + i +"/" + i + "_tab", dr=0.2)

#dossier_simulation = "../campagne_2Elements/G_c/4"
#liste_conf = sorted(os.listdir(dossier_simulation), key = num_key)
#liste_conf = [dossier_simulation + "/" + f for f in liste_conf if "conf" in f]
#save_concat_metrics(liste_conf, out_prefix="./Gc_4", dr=0.2)





#for file in liste_conf:
#    data = read_mfloe_conf(file)
#
#    # truc à extraire (???) à voir
    
#    res = compute_metrics(data, compute_gr=True, dr=0.2)
    







# HYPER IMPORTANT :
#colonne pour les floes : ['pos_x','pos_y','zpos','vel_x','vel_y','zvel','acc_x','acc_y','zacc', 'rot','vrot','arot','radius','height','inertia','mass']
#colonne pour les interactions :['i','j','isBonded','fn','fnb','ft','ftb','fs','fsb','A','coverage','Gc','dn0']


#data = read_mfloe_conf('../devel-manyelems/conf1')
#print("Meta:", data['meta'])
#print("Floe shape:", data['floe'].shape)          # (N,16)
#print("Interactions shape:", data['interactions'].shape)  # (M,13)
# accès à la colonne radius :
#radii = data['floe'][:, data['floe_cols'].index('radius')]
#print(radii)



