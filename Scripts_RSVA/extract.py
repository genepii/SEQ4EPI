import pandas as pd
import sys

def fix_malformed_csv(input_file, output_file, expected_columns):
    try:
        with open(input_file, 'r') as file:
            lines = file.readlines()

        # Vérifier chaque ligne et corriger les lignes mal formées
        fixed_lines = []
        for line in lines:
            # Compter le nombre de colonnes
            columns = line.split(',')
            if len(columns) == expected_columns:
                fixed_lines.append(line)
            else:
                print(f"Ligne mal formée trouvée: {line.strip()}")
                # Optionnel: tenter de corriger les lignes mal formées ici

        # Écrire les lignes corrigées dans un nouveau fichier
        with open(output_file, 'w') as file:
            file.writelines(fixed_lines)

        # Convertir le fichier corrigé avec des points-virgules comme séparateurs
        df = pd.read_csv(output_file, sep=',')
        df.to_csv(output_file, sep=';', index=False)
        
        print(f"Fichier converti et sauvegardé dans {output_file}")
    except Exception as e:
        print(f"Erreur pendant la conversion: {e}")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python fix_malformed_csv.py <input_file> <output_file> <expected_columns>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    expected_columns = int(sys.argv[3])

    fix_malformed_csv(input_file, output_file, expected_columns)



