# genui/src/genui/generators/extensions/genuireinvent/genuimodels/algorithms.py

from __future__ import annotations

from abc import ABC
import pickle

from genui.models.genuimodels import bases
from genui.models.models import ModelFileFormat


class ReinventAlgorithm(bases.Algorithm, ABC):
    """
    Thin adapter around the external REINVENT CLI.

    Training is delegated to the model facade returned by ReinventNet.getModel().
    The builder/GenUI training pipeline will call Algorithm.fit(); internally that
    triggers the CLI-based transfer learning in ReinventNet.run_transfer_learning().
    """

    def __init__(self, builder, callback=None):
        super().__init__(builder, callback)
        self.train_params = {}
        self._model = None

    @classmethod
    def getFileFormats(cls, attach_to=None):
        pkg, _ = ModelFileFormat.objects.get_or_create(
            fileExtension=".pkg",
            defaults={
                "description": "Serialized metadata for REINVENT (e.g., produced checkpoint path).",
            },
        )
        if attach_to:
            cls.attachToInstance(attach_to, [pkg], attach_to.fileFormats)

    @classmethod
    def getModes(cls):
        return [cls.GENERATOR]

    @property
    def model(self):
        return self._model

    def predict(self, X):
        # REINVENT net here is used as a generator; prediction isn't applicable.
        return [], None

    def sample(self, n_samples, from_inputs=None):
        raise NotImplementedError("Sampling is not implemented for the REINVENT CLI adapter.")

    def getSerializer(self):
        """
        Persist minimal metadata needed to re-associate a trained artifact with this Algorithm.

        NOTE:
        - The real checkpoint is managed by ReinventNet.checkpointFile (AUX ModelFile).
        - This .pkg payload is primarily for GenUI's generic model-serialization contract.
        """
        def _save(path: str):
            payload = {}
            try:
                # facade.getModel() returns {"checkpoint": "..."} in your setup
                payload = self.model.getModel() if self.model else {}
            except Exception:
                payload = {}
            with open(path, "wb") as f:
                pickle.dump(payload, f, protocol=pickle.HIGHEST_PROTOCOL)
        return _save

    def getDeserializer(self):
        """
        Restore the facade. Nothing is loaded into RAM.
        """
        def _load(path: str):
            try:
                with open(path, "rb") as f:
                    _ = pickle.load(f)  # kept for compatibility; optional
            except Exception:
                pass

            # facade may be lazily created by subclasses; if present, let it no-op load.
            if self.model:
                try:
                    self.model.loadStatesFromFile(path)
                except Exception:
                    pass
            return self.model
        return _load


class ReinventNetwork(ReinventAlgorithm):
    name = "ReinventNet"

    def __init__(self, builder, callback=None):
        super().__init__(builder, callback)
        # builder.instance is a ReinventNet (Django model) which returns the CLI facade
        self._model = self.builder.instance.getModel()

    def fit(self, X=None, y=None):
        # Refresh facade (safe in case builder.instance mutated)
        self._model = self.builder.instance.getModel()

        # Triggers: prepareData() + run_transfer_learning() (subprocess)
        self._model.fit(X=X, y=y)

        # Inform pipeline that "an epoch-like thing happened"
        if self.callback:
            self.callback(None)

        return self


class ReinventSampler:
    """
    Utility class for sampling SMILES from a trained REINVENT model.
    Used by Reinvent.get() to generate molecules on-demand.
    """

    def __init__(self, model_path, device="cpu"):
        """
        Initialize the sampler with a trained model.

        :param model_path: Path to the .ckpt or .prior model file
        :param device: Device to run on ('cpu' or 'cuda')
        """
        import torch

        self.device = device
        self.model_path = model_path

        # Load the model using torch.load (works for REINVENT checkpoints)
        try:
            # Load checkpoint
            # Note: weights_only=False is required for REINVENT models (PyTorch 2.6+)
            checkpoint = torch.load(model_path, map_location=device, weights_only=False)

            # Extract the model (checkpoint format may vary)
            if isinstance(checkpoint, dict):
                if 'model' in checkpoint:
                    self.model = checkpoint['model']
                elif 'network' in checkpoint:
                    self.model = checkpoint['network']
                elif 'model_state_dict' in checkpoint:
                    # Need to reconstruct model from state dict
                    # For now, just use the whole checkpoint
                    self.model = checkpoint
                else:
                    # Assume the checkpoint itself is the model
                    self.model = checkpoint
            else:
                # Checkpoint is the model directly
                self.model = checkpoint

            # Move to device if it's a PyTorch module
            if hasattr(self.model, 'to'):
                self.model = self.model.to(device)
                self.model.eval()

        except Exception as e:
            raise RuntimeError(f"Could not load REINVENT model from {model_path}: {str(e)}")

    def sample(self, n_samples):
        """
        Sample SMILES from the model.

        :param n_samples: Number of SMILES to generate
        :return: List of SMILES strings
        """
        import torch

        smiles_list = []
        batch_size = min(128, n_samples)

        with torch.no_grad():
            while len(smiles_list) < n_samples:
                current_batch_size = min(batch_size, n_samples - len(smiles_list))

                try:
                    # Try different sampling methods
                    if hasattr(self.model, 'sample_smiles'):
                        # REINVENT API with direct SMILES output
                        batch_smiles = self.model.sample_smiles(current_batch_size)
                    elif hasattr(self.model, 'sample'):
                        # Standard sample method - returns (seqs, smiles, nlls)
                        result = self.model.sample(current_batch_size)
                        if isinstance(result, tuple) and len(result) >= 2:
                            batch_smiles = result[1]  # SMILES are second element
                        else:
                            batch_smiles = result
                    else:
                        raise NotImplementedError(
                            f"Model does not have a sample method. Available methods: {dir(self.model)}"
                        )

                    # Convert to list if needed
                    if isinstance(batch_smiles, str):
                        batch_smiles = [batch_smiles]
                    elif not isinstance(batch_smiles, list):
                        # Try to convert to list
                        batch_smiles = list(batch_smiles)

                    smiles_list.extend(batch_smiles[:current_batch_size])

                except Exception as e:
                    raise RuntimeError(f"Error during sampling: {str(e)}")

        return smiles_list[:n_samples]

