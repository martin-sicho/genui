"""
helpers

Created by: Martin Sicho
On: 02-03-20, 14:21
"""
from io import StringIO, BytesIO

from django.core.files.base import ContentFile
from django.db import transaction

from compounds.models import PictureFormat, Molecule, MoleculePic
from rdkit.Chem import Draw, rdDepictor
from rdkit import Chem


def createPic(mol : Molecule, format : PictureFormat):
    content = None
    rd_mol = mol.molObject
    if format.extension == '.svg':
        rdDepictor.Compute2DCoords(rd_mol)
        mc = Chem.Mol(rd_mol.ToBinary())
        Chem.Kekulize(mc)
        drawer = Draw.MolDraw2DSVG(400,400)
        drawer.DrawMolecule(mc)
        drawer.FinishDrawing()
        content = StringIO()
        content.write(drawer.GetDrawingText())
    elif format.extension == '.png':
        rdDepictor.Compute2DCoords(rd_mol)
        content = Draw.MolToImage(rd_mol, size=(400, 400))
        content_io = BytesIO()
        content.save(content_io, format='PNG')
        content = ContentFile(content_io.getvalue())

    with transaction.atomic():
        img = MoleculePic.objects.create(
                molecule=mol,
                format=format
            )
        img.image.save(f'mol_{mol.id}{format.extension}', content)
        img.save()