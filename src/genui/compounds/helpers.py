"""
helpers

Created by: Martin Sicho
On: 02-03-20, 14:21
"""
from io import StringIO, BytesIO

from django.core.files.base import ContentFile
from django.db import transaction

from . import models
from rdkit.Chem import Draw, rdDepictor
from rdkit import Chem


def createPic(mol : "models.Molecule", format : "models.PictureFormat"):
    content = None
    rd_mol = mol.rdMol
    if format.extension == '.svg':
        rdDepictor.Compute2DCoords(rd_mol)
        mc = Chem.Mol(rd_mol.ToBinary())
        Chem.Kekulize(mc)
        drawer = Draw.MolDraw2DSVG(400,400)
        drawer.DrawMolecule(mc)
        drawer.FinishDrawing()
        content_io = StringIO()
        content_io.write(drawer.GetDrawingText())
    elif format.extension == '.png':
        rdDepictor.Compute2DCoords(rd_mol)
        content = Draw.MolToImage(rd_mol, size=(400, 400))
        content_io = BytesIO()
        content.save(content_io, format='PNG')
    else:
        raise NotImplementedError('Unsupported file format:', format.extension)

    content = ContentFile(content_io.getvalue())
    with transaction.atomic():
        img = models.MoleculePic.objects.create(
                molecule=mol,
                format=format
            )
        img.image.save(f'mol_{mol.id}{format.extension}', content)
        img.save()

    return img