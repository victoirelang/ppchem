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




