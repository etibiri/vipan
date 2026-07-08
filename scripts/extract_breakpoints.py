import json
import csv
import sys
import os

def extract_gard_results(json_path, csv_path):
    if not os.path.exists(json_path):
        print(f"Erreur : {json_path} non trouvé.")
        sys.exit(1)

    with open(json_path, 'r') as f:
        data = json.load(f)

    # Extraction des points de cassure (breakpoints)
    # GARD stocke les modèles testés, on cherche celui avec le meilleur AIC
    bps = data.get("breakpoints", [])
    
    with open(csv_path, 'w', newline='') as f_out:
        writer = csv.writer(f_out)
        writer.writerow(["Breakpoint_Index", "Nucleotide_Position", "AIC_Score"])
        
        if not bps:
            print("Aucun breakpoint détecté par GARD.")
            return

        # On récupère les positions du meilleur modèle
        # Souvent le dernier ou celui avec le plus de breakpoints validés
        last_model = bps[-1] 
        positions = last_model.get("breakpoints", [])
        aic = last_model.get("AIC", "N/A")

        for i, pos in enumerate(positions):
            writer.writerow([i + 1, pos, aic])

    print(f"Succès : {len(positions)} breakpoints extraits dans {csv_path}")

if __name__ == "__main__":
    extract_gard_results(sys.argv[1], sys.argv[2])
