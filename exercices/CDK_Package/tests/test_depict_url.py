import urllib.parse
import requests
from IPython.display import SVG, display

cdk_depict_link = 'https://www.simolecule.com/cdkdepict/depict/bow'


def smiles_depict_url(smiles: str, format: str = 'svg') -> str:
    """
    Generate the URL for the depiction of a SMILES string.
    Args:
        smiles: smiles string to depict
        format: 'svg', 'pdf', 'png', etc.
        use_internal_service: whether to use the service deployed on ZC2 (True)
            or the one available on the simolecule website (False).
    Returns:
        URL string
    """
    rendering_service = cdk_depict_link
    params = {
        'smi': smiles,
        'zoom': '1.0',
        'abbr': 'on',
        'hdisp': 'bridgehead',
        'showtitle': 'false',
        'annotate': 'none'
    }
    params_str = urllib.parse.urlencode(params)
    return f'{rendering_service}/{format}?{params_str}'

def display_svg(url: str) -> None:
    response = requests.get(url)
    if response.status_code == 200:
        svg_data = response.text
        display(SVG(svg_data))
    else:
        print("Failed to retrieve SVG: Status code", response.status_code)

smiles = 'CCOC(=O)C1=CC=CC=C1C(=O)OCC' 
url = smiles_depict_url(smiles)
display_svg(url)