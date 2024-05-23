def visualize_molecules_for_cream3(df, cream_name):
    """
    Affiche les molécules d'une certaine catégorie de crèmes avec des atomes spécifiques mis en évidence.
    
    Parameters:
    df (pd.DataFrame): DataFrame contenant les données.
    cream_name (str): Nom de la crème à visualiser.
    """
    # Filtrer les molécules pour la catégorie de crème spécifiée
    filtered_df = df[df['Cream'] == cream_name]
    # Supprimer les lignes avec des valeurs manquantes dans 'Smiles'
    filtered_df = filtered_df.dropna(subset=['Smiles'])
    # Supprimer les lignes avec des valeurs non-string dans 'Smiles'
    filtered_df = filtered_df[filtered_df['Smiles'].apply(lambda x: isinstance(x, str))]
    
    # Définir les couleurs pour les atomes spécifiques
    color_map = {
        'Isopropyl': (0, 1, 0),  # Vert
        'NaOH': (0, 0, 1),       # Bleu
        'Linalool': (0.6, 0, 0.6),   # Violet 1
        'Citronellol': (0.7, 0, 0.7), # Violet 2
        'Limonene': (0.8, 0, 0.8),   # Violet 3
        'Benzyl Alcohol': (0.9, 0, 0.9), # Violet 4
        'Benzyl Salicylate': (1, 0, 1)   # Violet 5
    }
    
    # Définir les motifs SMARTS pour les composés spécifiques
    smarts_patterns = {
        'Isopropyl': 'CC(C)O',
        'NaOH': '[Na+].[OH-]',
        'Linalool': 'CC(=C)CCC1=CC=C(C=C1)O',
        'Citronellol': 'CC(C)CCC(C)CO',
        'Limonene': 'CC1=CC2CCC1(C)C=C2',
        'Benzyl Alcohol': 'C1=CC=C(C=C1)CO',
        'Benzyl Salicylate': 'C1=CC=C(C=C1)COC(=O)C2=CC=CC=C2'
    }
    
    mols = []
    atom_colors = []
    
    # Parcourir les molécules filtrées et préparer les images
    for index, row in filtered_df.iterrows():
        smi = row['Smiles']
        mol = Chem.MolFromSmiles(smi)
        
        # Ignorer les SMILES invalides
        if mol is None:
            continue
        
        highlight_dict = {}
        # Identifier les atomes à mettre en évidence
        for compound, smarts in smarts_patterns.items():
            pattern = Chem.MolFromSmarts(smarts)
            matches = mol.GetSubstructMatches(pattern)
            for match in matches:
                for idx in match:
                    highlight_dict[idx] = color_map[compound]
        
        # Ajouter la molécule et le dictionnaire de surbrillance à la liste
        mols.append(mol)
        atom_colors.append(highlight_dict)
    
    # Afficher chaque molécule individuellement
    for mol, highlights in zip(mols, atom_colors):
        drawer = rdMolDraw2D.MolDraw2DSVG(500, 500)
        drawer.drawOptions().useBWAtomPalette()
        drawer.DrawMolecule(mol, highlightAtoms=list(highlights.keys()), highlightAtomColors=highlights)
        drawer.FinishDrawing()
        
        # Afficher l'image de la molécule
        svg = drawer.GetDrawingText().replace('svg:', '')
        display(SVG(svg))














def check_smiles_validity(df, smiles_column='Smiles'):
    """
    Vérifie la validité des chaînes SMILES dans un DataFrame.
    
    Parameters:
    df (pd.DataFrame): DataFrame contenant les données.
    smiles_column (str): Nom de la colonne contenant les chaînes SMILES.
    
    Returns:
    dict: Dictionnaire avec les résultats de la vérification. Les clés sont les indices des lignes du DataFrame,
          et les valeurs sont des booléens indiquant si la chaîne SMILES est valide (True) ou non (False).
    """
    validity_dict = {}
    
    for index, row in df.iterrows():
        smi = row[smiles_column]
        if not isinstance(smi, str):
            validity_dict[index] = False
            continue
        
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            validity_dict[index] = False
        else:
            validity_dict[index] = True
    
    return validity_dict















def visualize_molecules_for_cream1(df, cream_name):
    """
    Affiche les molécules d'une certaine catégorie de crèmes avec des atomes spécifiques mis en évidence.
    
    Parameters:
    df (pd.DataFrame): DataFrame contenant les données.
    cream_name (str): Nom de la crème à visualiser.
    """
    # Filtrer les molécules pour la catégorie de crème spécifiée
    filtered_df = df[df['Cream'] == cream_name]
    # Supprimer les lignes avec des valeurs manquantes dans 'Smiles'
    filtered_df = filtered_df.dropna(subset=['Smiles'])
    # Supprimer les lignes avec des valeurs non-string dans 'Smiles'
    filtered_df = filtered_df[filtered_df['Smiles'].apply(lambda x: isinstance(x, str))]
    
    # Définir les couleurs pour les atomes spécifiques
    color_map = {
        'Isopropyl': (0, 1, 0),  # Vert
        'NaOH': (0, 0, 1),       # Bleu
        'Linalool': (0.6, 0, 0.6),   # Violet 1
        'Citronellol': (0.7, 0, 0.7), # Violet 2
        'Limonene': (0.8, 0, 0.8),   # Violet 3
        'Benzyl Alcohol': (0.9, 0, 0.9), # Violet 4
        'Benzyl Salicylate': (1, 0, 1)   # Violet 5
    }
    
    # Définir les motifs SMARTS pour les composés spécifiques
    smarts_patterns = {
        'Isopropyl': 'CC(C)O',
        'NaOH': '[Na+].[OH-]',
        'Linalool': 'CC(C)CC1=CC=C(C=C1)O',
        'Citronellol': 'CC(C)CCC(C)(C)CO',
        'Limonene': 'CC1=CC=CCC1(C)C',
        'Benzyl Alcohol': 'C1=CC=C(C=C1)CO',
        'Benzyl Salicylate': 'C1=CC=C(C=C1)OC(=O)C2=CC=CC=C2'
    }
    
    mols = []
    atom_colors = []
    
    # Parcourir les molécules filtrées et préparer les images
    for index, row in filtered_df.iterrows():
        smi = row['Smiles']
        try:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                raise ValueError("Mol is None")
        except Exception as e:
            print(f"Erreur de parsing SMILES pour: {smi}, erreur: {e}")
            continue

        highlight_dict = {}
        # Identifier les atomes à mettre en évidence
        for compound, smarts in smarts_patterns.items():
            pattern = Chem.MolFromSmarts(smarts)
            matches = mol.GetSubstructMatches(pattern)
            for match in matches:
                for idx in match:
                    highlight_dict[idx] = color_map[compound]
        
        # Ajouter la molécule et le dictionnaire de surbrillance à la liste
        mols.append(mol)
        atom_colors.append(highlight_dict)
    
    # Générer la grille d'images sans légendes et avec des images plus grandes
    n_mols = len(mols)
    n_cols = 3
    n_rows = (n_mols + n_cols - 1) // n_cols  # Calculer le nombre de lignes nécessaire
    mol_size = (500, 500)
    
    drawer = rdMolDraw2D.MolDraw2DSVG(n_cols * mol_size[0], n_rows * mol_size[1])
    drawer.SetFontSize(1.0)
    
    for i, mol in enumerate(mols):
        row, col = divmod(i, n_cols)
        drawer.SetOffset(col * mol_size[0], row * mol_size[1])
        
        # Dessiner chaque molécule individuellement
        drawer.DrawMolecule(mol, highlightAtoms=list(atom_colors[i].keys()), highlightAtomColors=atom_colors[i])
    
    drawer.FinishDrawing()
    
    # Afficher l'image de la grille
    svg = drawer.GetDrawingText().replace('svg:', '')
    display(SVG(svg))




def visualize_molecules_for_cream4(df, cream_name):
    """
    Affiche les molécules d'une certaine catégorie de crèmes avec des atomes spécifiques mis en évidence.
    
    Parameters:
    df (pd.DataFrame): DataFrame contenant les données.
    cream_name (str): Nom de la crème à visualiser.
    """
    # Filtrer les molécules pour la catégorie de crème spécifiée
    filtered_df = df[df['Cream'] == cream_name]
    # Supprimer les lignes avec des valeurs manquantes dans 'Smiles'
    filtered_df = filtered_df.dropna(subset=['Smiles'])
    # Supprimer les lignes avec des valeurs non-string dans 'Smiles'
    filtered_df = filtered_df[filtered_df['Smiles'].apply(lambda x: isinstance(x, str))]
    
    # Définir les couleurs pour les atomes spécifiques
    color_map = {
        'Isopropyl': (0, 1, 0),  # Vert
        'NaOH': (0, 0, 1),       # Bleu
        'Linalool': (0.6, 0, 0.6),   # Violet 1
        'Citronellol': (0.7, 0, 0.7), # Violet 2
        'Limonene': (0.8, 0, 0.8),   # Violet 3
        'Benzyl Alcohol': (0.9, 0, 0.9), # Violet 4
        'Benzyl Salicylate': (1, 0, 1)   # Violet 5
    }
    
    # Définir les motifs SMARTS pour les composés spécifiques
    smarts_patterns = {
        'Isopropyl': 'CC(C)O',
        'NaOH': '[Na+].[OH-]',
        'Linalool': 'CC(=C)CCC1=CC=C(C=C1)O',
        'Citronellol': 'CC(C)CCC(C)CO',
        'Limonene': 'CC1=CC2CCC1(C)C=C2',
        'Benzyl Alcohol': 'C1=CC=C(C=C1)CO',
        'Benzyl Salicylate': 'C1=CC=C(C=C1)COC(=O)C2=CC=CC=C2'
    }
    
    mols = []
    atom_colors = []
    
    # Parcourir les molécules filtrées et préparer les images
    for index, row in filtered_df.iterrows():
        smi = row['Smiles']
        mol = Chem.MolFromSmiles(smi)
        
        # Ignorer les SMILES invalides
        if mol is None:
            continue
        
        highlight_dict = {}
        # Identifier les atomes à mettre en évidence
        for compound, smarts in smarts_patterns.items():
            pattern = Chem.MolFromSmarts(smarts)
            matches = mol.GetSubstructMatches(pattern)
            for match in matches:
                for idx in match:
                    highlight_dict[idx] = color_map[compound]
        
        # Ajouter la molécule et le dictionnaire de surbrillance à la liste
        mols.append(mol)
        atom_colors.append(highlight_dict)
    
    # Afficher chaque molécule individuellement
    for mol, highlights in zip(mols, atom_colors):
        drawer = rdMolDraw2D.MolDraw2DSVG(500, 500)
        drawer.drawOptions().useBWAtomPalette()
        drawer.DrawMolecule(mol, highlightAtoms=list(highlights.keys()), highlightAtomColors=highlights)
        drawer.FinishDrawing()
        
        # Afficher l'image de la molécule
        svg = drawer.GetDrawingText().replace('svg:', '')
        display(SVG(svg))


#'Linalool': 'CC(=C)CCC1=CC=C(C=C1)O',
       # 'Citronellol': 'CC(C)CCC(C)CO',
       # 'Limonene': 'CC1=CC2CCC1(C)C=C2',
       # 'Benzyl Alcohol': 'C1=CC=C(C=C1)CO',
       # 'Benzyl Salicylate': 'C1=CC=C(C=C1)COC(=O)C2=CC=CC=C2'






import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D
from IPython.display import SVG, display

def visualize_molecules_for_cream_smarts(df, cream_name):
    """
    Affiche les molécules d'une certaine catégorie de crèmes avec des atomes spécifiques mis en évidence.
    
    Parameters:
    df (pd.DataFrame): DataFrame contenant les données.
    cream_name (str): Nom de la crème à visualiser.
    """
    # Filtrer les molécules pour la catégorie de crème spécifiée
    filtered_df = df[df['Cream'] == cream_name]
    # Supprimer les lignes avec des valeurs manquantes dans 'Smiles'
    filtered_df = filtered_df.dropna(subset=['Smiles'])
    # Supprimer les lignes avec des valeurs non-string dans 'Smiles'
    filtered_df = filtered_df[filtered_df['Smiles'].apply(lambda x: isinstance(x, str))]
    
    print(f"Total molécules pour la crème {cream_name}: {len(filtered_df)}")
    
    # Définir les couleurs pour les atomes spécifiques
    color_map = {
        'Isopropyl': (0, 1, 0),  # Vert
        'NaOH': (0, 0, 1),       # Bleu
        'Linalool': (0.25, 0, 0.25),   # Violet 1
        'Citronellol': (0.4, 0, 0.4), # Violet 2
        'Limonene': (0.6, 0, 0.6),   # Violet 3
        'Benzyl Alcohol': (0.8, 0, 0.8), # Violet 4
        'Benzyl Salicylate': (1, 0, 1)   # Violet 5
    }
    
    # Définir les motifs SMARTS pour les composés spécifiques
    smarts_patterns = {
        'Isopropyl': '[CH3][CH](O)[CH3]',
        'NaOH': '[Na+].[OH-]',
        'Linalool': 'CC(=C)CCCc1ccccc1O',
        'Citronellol': 'CC(C)CCC(C)CO',
        'Limonene': 'CC1=CC2CCC1(C)C=C2',
        'Benzyl Alcohol': 'c1ccccc1CO',  
        'Benzyl Salicylate': 'c1ccc(cc1)COC(=O)c2ccccc2'
    }
    
    mols = []
    atom_colors = []
    
    # Parcourir les molécules filtrées et préparer les images
    for index, row in filtered_df.iterrows():
        smi = row['Smiles']
        mol = Chem.MolFromSmiles(smi)
        
        # Ignorer les SMILES invalides
        if mol is None:
            continue
        
        highlight_dict = {}
        # Identifier les atomes à mettre en évidence
        for compound, smarts in smarts_patterns.items():
            pattern = Chem.MolFromSmarts(smarts)
            matches = mol.GetSubstructMatches(pattern)
            if matches:
                print(f"Match trouvé pour {compound} dans la molécule {smi}")
            for match in matches:
                for idx in match:
                    highlight_dict[idx] = color_map[compound]
        
        # Ajouter la molécule et le dictionnaire de surbrillance à la liste
        mols.append(mol)
        atom_colors.append(highlight_dict)
    
    # Afficher chaque molécule individuellement
    for mol, highlights in zip(mols, atom_colors):
        drawer = rdMolDraw2D.MolDraw2DSVG(500, 500)
        drawer.drawOptions().useBWAtomPalette()
        drawer.DrawMolecule(mol, highlightAtoms=list(highlights.keys()), highlightAtomColors=highlights)
        drawer.FinishDrawing()
        
        # Afficher l'image de la molécule
        svg = drawer.GetDrawingText().replace('svg:', '')
        display(SVG(svg))

















def visualize_molecules_for_cream1(df, cream_name):
    """
    Affiche les molécules d'une certaine catégorie de crèmes avec des atomes spécifiques mis en évidence.
    
    Parameters:
    df (pd.DataFrame): DataFrame contenant les données.
    cream_name (str): Nom de la crème à visualiser.
    """
    # Filtrer les molécules pour la catégorie de crème spécifiée
    filtered_df = df[df['Cream'] == cream_name]
    # Supprimer les lignes avec des valeurs manquantes dans 'Smiles'
    filtered_df = filtered_df.dropna(subset=['Smiles'])
    # Supprimer les lignes avec des valeurs non-string dans 'Smiles'
    filtered_df = filtered_df[filtered_df['Smiles'].apply(lambda x: isinstance(x, str))]
    
    print(f"Total molécules pour la crème {cream_name}: {len(filtered_df)}")
    
    # Définir les couleurs pour les atomes spécifiques
    color_map = {
        'Isopropyl': (0, 1, 0),  # Vert
        'NaOH': (0, 0, 1),       # Bleu
        'Linalool': (0.6, 0, 0.6),   # Violet 1
        'Citronellol': (0.7, 0, 0.7), # Violet 2
        'Limonene': (0.8, 0, 0.8),   # Violet 3
        'Benzyl Alcohol': (0.9, 0, 0.9), # Violet 4
        'Benzyl Salicylate': (1, 0, 1)   # Violet 5
    }
    
    # Définir les motifs SMARTS pour les composés spécifiques
    smarts_patterns = {
        'Isopropyl': '[CH3][CH](O)[CH3]',
        'NaOH': '[Na+].[OH-]',
        'Linalool': 'CC(=C)CCC1=CC=C(C=C1)O',
        'Citronellol': 'CC(C)CCC(C)CO',
        'Limonene': 'CC1=CC[C@H](C=C1)C(C)C',
        'Benzyl Alcohol': 'C1=CC=C(C=C1)CO',
        'Benzyl Salicylate': 'COC1=CC=CC=C1C(=O)O'
    }
    
    mols = []
    atom_colors = []
    
    # Parcourir les molécules filtrées et préparer les images
    for index, row in filtered_df.iterrows():
        smi = row['Smiles']
        mol = Chem.MolFromSmiles(smi)
        
        # Ignorer les SMILES invalides
        if mol is None:
            print(f"SMILES invalide: {smi}")
            continue
        
        highlight_dict = {}
        # Identifier les atomes à mettre en évidence
        for compound, smarts in smarts_patterns.items():
            pattern = Chem.MolFromSmarts(smarts)
            matches = mol.GetSubstructMatches(pattern)
            if matches:
                print(f"Match trouvé pour {compound} dans la molécule {smi}")
            for match in matches:
                for idx in match:
                    highlight_dict[idx] = color_map[compound]
        
        # Ajouter la molécule et le dictionnaire de surbrillance à la liste
        mols.append(mol)
        atom_colors.append(highlight_dict)
    
    # Afficher chaque molécule individuellement
    for mol, highlights in zip(mols, atom_colors):
        drawer = rdMolDraw2D.MolDraw2DSVG(500, 500)
        drawer.drawOptions().useBWAtomPalette()
        drawer.DrawMolecule(mol, highlightAtoms=list(highlights.keys()), highlightAtomColors=highlights)
        drawer.FinishDrawing()
        
        # Afficher l'image de la molécule
        svg = drawer.GetDrawingText().replace('svg:', '')
        display(SVG(svg))




def visualize_molecules_for_cream(df, cream_name):
    """
    Affiche les molécules d'une certaine catégorie de crèmes avec des atomes spécifiques mis en évidence.
    
    Parameters:
    df (pd.DataFrame): DataFrame contenant les données.
    cream_name (str): Nom de la crème à visualiser.
    """
    # Filtrer les molécules pour la catégorie de crème spécifiée
    filtered_df = df[df['Cream'] == cream_name]
    # Supprimer les lignes avec des valeurs manquantes dans 'Smiles'
    filtered_df = filtered_df.dropna(subset=['Smiles'])
    # Supprimer les lignes avec des valeurs non-string dans 'Smiles'
    filtered_df = filtered_df[filtered_df['Smiles'].apply(lambda x: isinstance(x, str))]
    
    # Définir les couleurs pour les atomes spécifiques
    color_map = {
        'Isopropyl': (0, 1, 0),  # Vert
        'NaOH': (0, 0, 1),       # Bleu
        'Linalool': (0.25, 0, 0.25),   # Violet 1
        'Citronellol': (0.4, 0, 0.4), # Violet 2
        'Limonene': (0.6, 0, 0.6),   # Violet 3
        'Benzyl Alcohol': (0.8, 0, 0.8), # Violet 4
        'Benzyl Salicylate': (1, 0, 1)   # Violet 5
    }
    
    # Définir les motifs SMARTS pour les composés spécifiques
    smarts_patterns = {
        'Isopropyl': '[CH3][CH](O)[CH3]',
        'NaOH': '[Na+].[OH-]',
        'Linalool': 'CC(=C)CCC1=CC=C(C=C1)O',
        'Citronellol': 'CC(C)CCC(C)CO',
        'Limonene': 'CC1=CC2CCC1(C)C=C2',
        'Benzyl Alcohol': 'C1=CC=C(C=C1)CO',  
        'Benzyl Salicylate': 'C1=CC=C(C=C1)COC(=O)C2=CC=CC=C2'
    }
    
    mols = []
    atom_colors = []
    
    # Parcourir les molécules filtrées et préparer les images
    for index, row in filtered_df.iterrows():
        smi = row['Smiles']
        mol = Chem.MolFromSmiles(smi)
        
        # Ignorer les SMILES invalides
        if mol is None:
            continue
        
        highlight_dict = {}
        # Identifier les atomes à mettre en évidence
        for compound, smarts in smarts_patterns.items():
            pattern = Chem.MolFromSmarts(smarts)
            matches = mol.GetSubstructMatches(pattern)
            if matches:
                print(f"Match trouvé pour {compound} dans la molécule {smi}")
            for match in matches:
                for idx in match:
                    highlight_dict[idx] = color_map[compound]
        
        # Ajouter la molécule et le dictionnaire de surbrillance à la liste
        mols.append(mol)
        atom_colors.append(highlight_dict)
    
    # Afficher chaque molécule individuellement
    for mol, highlights in zip(mols, atom_colors):
        drawer = rdMolDraw2D.MolDraw2DSVG(500, 500)
        drawer.drawOptions().useBWAtomPalette()
        drawer.DrawMolecule(mol, highlightAtoms=list(highlights.keys()), highlightAtomColors=highlights)
        drawer.FinishDrawing()
        
        # Afficher l'image de la molécule
        svg = drawer.GetDrawingText().replace('svg:', '')
        display(SVG(svg))

