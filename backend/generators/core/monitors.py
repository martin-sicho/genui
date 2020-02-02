"""
monitors

Created by: Martin Sicho
On: 30-01-20, 15:05
"""
import numpy as np

from drugex.api.agent.callbacks import AgentMonitor
from drugex.api.model.callbacks import PretrainingMonitor
from generators import models
from generators.core import metrics


class DrugExNetMonitor(PretrainingMonitor):

    def __init__(self, builder, original_callback=None):
        self.original_call = original_callback
        self.builder = builder
        self.loss_train = None
        self.loss_valid = None
        self.error_rate = None
        self.best_error = None
        self.current_step = None
        self.current_epoch = None
        self.total_steps = None
        self.total_epochs = None
        self.best_state = None
        self.best_yet = False
        self.current_model = None

    def savePerformance(self, metric, value, isValidation, note=""):
        return models.ModelPerformanceDrugEx.objects.create(
            metric=metric,
            value=value,
            isOnValidationSet=isValidation,
            model=self.builder.instance,
            epoch=self.current_epoch,
            step=self.current_step,
            note=note
        )

    @property
    def last_step(self):
        # TODO: determine where the training left off last time

        return 0

    @property
    def last_epoch(self):
        # TODO: determine where the training left off last time

        return 0

    def finalizeStep(self, current_epoch: int, current_batch: int, current_step: int, total_epochs: int,
                     total_batches: int, total_steps: int):

        self.current_step = current_step + self.last_step
        self.current_epoch = current_epoch + self.last_epoch
        self.total_steps = total_steps + self.last_step
        self.total_epochs = total_epochs + self.last_epoch

        info = "Epoch: %d step: %d error_rate: %.3f loss_train: %.3f loss_valid %.3f" % (self.current_epoch, self.current_step, self.error_rate, self.loss_train, self.loss_valid if self.loss_valid is not None else np.inf)
        print(info)

        if self.loss_train is not None:
            self.savePerformance(metrics.DrugExLoss.getDjangoModel(), self.loss_train, False)
        if self.loss_valid is not None:
            self.savePerformance(metrics.DrugExLoss.getDjangoModel(), self.loss_valid, True)
        if self.error_rate:
            if self.best_yet:
                self.savePerformance(metrics.SMILESErrorRate.getDjangoModel(), self.error_rate, False, note=f"Minimum {metrics.SMILESErrorRate.name}, yet.")
            else:
                self.savePerformance(metrics.SMILESErrorRate.getDjangoModel(), self.error_rate, False)

        if self.original_call:
            self.original_call(self)

    def performance(self, loss_train, loss_valid, error_rate, best_error):
        self.loss_train = loss_train
        self.loss_valid = loss_valid
        self.error_rate = error_rate
        self.best_error = best_error

    def smiles(self, smiles, is_valid):
        pass

    def model(self, model):
        self.current_model = model

    def state(self, current_state, is_best=False):
        self.best_yet = is_best
        if self.best_yet:
            self.best_state = current_state
            self.builder.saveFile()

    def close(self):
        pass

    def getState(self):
        return self.best_state

class DrugExAgentMonitor(AgentMonitor):

    def __init__(self, builder, original_callback=None):
        self.builder = builder
        self.original_callback = original_callback
        self.best_yet = None
        self.best_state = None

    def finalizeEpoch(self, current_epoch, total_epochs):
        pass

        if self.original_callback:
            self.original_callback(self)

    def smiles(self, smiles, score):
        pass

    def performance(self, scores, valids, criterion, best_score):
        pass

    def model(self, model):
        pass

    def state(self, current_state, is_best=False):
        self.best_yet = is_best
        if self.best_yet:
            self.best_state = current_state
            self.builder.saveFile()

    def close(self):
        pass

    def getState(self):
        return self.best_state