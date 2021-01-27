"""
monitors

Created by: Martin Sicho
On: 30-01-20, 15:05
"""
import weakref

import numpy as np

from drugex.api.agent.callbacks import AgentMonitor
from drugex.api.model.callbacks import PretrainingMonitor
from . import metrics
from ..models import ModelPerformanceDrugEx, ModelPerformanceDrugExAgent


class DrugExNetMonitor(PretrainingMonitor):

    def __init__(self, builder, original_callback=None):
        self.original_call = original_callback
        self._builder = weakref.ref(builder)
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

    @property
    def builder(self):
        ret = self._builder()
        if ret:
            return ret
        else:
            raise LookupError("Builder was destroyed before being referenced!")

    def savePerformance(self, metric, value, isValidation, note=""):
        return ModelPerformanceDrugEx.objects.create(
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

        if current_epoch != self.current_epoch:
            self.builder.recordProgress()

        self.current_step = current_step
        self.current_epoch = current_epoch
        self.total_steps = total_steps + self.last_step
        self.total_epochs = total_epochs + self.last_epoch

        info = "Epoch: %d step: %d error_rate: %.3f loss_train: %.3f loss_valid %.3f" % (self.current_epoch, self.current_step, self.error_rate, self.loss_train, self.loss_valid if self.loss_valid is not None else np.inf)
        print(info)

        if self.loss_train is not None:
            self.savePerformance(metrics.DrugExLoss.getDjangoModel(), self.loss_train, False)
        if self.loss_valid is not None:
            self.savePerformance(metrics.DrugExLoss.getDjangoModel(), self.loss_valid, True)
        if self.error_rate is not None:
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
        self._builder = weakref.ref(builder)
        self.original_callback = original_callback
        self.best_yet = None
        self.best_state = None
        self.current_epoch = None
        self.current_criterion = None
        self.current_error = None
        self.current_mean_score = None

    @property
    def builder(self):
        ret = self._builder()
        if ret:
            return ret
        else:
            raise LookupError("Builder was destroyed before being referenced!")

    @property
    def last_epoch(self):
        # TODO: determine where the training left off last time

        return 0

    def savePerformance(self, metric, value, note=""):
        return ModelPerformanceDrugExAgent.objects.create(
            metric=metric,
            value=value,
            model=self.builder.instance,
            epoch=self.current_epoch,
            note=note
        )

    def finalizeEpoch(self, current_epoch, total_epochs):
        self.current_epoch = current_epoch + self.last_epoch + 1

        print("Epoch+: %d average: %.4f valid: %.4f unique: %.4f" % (self.current_epoch, self.current_mean_score, self.current_error, self.current_criterion))

        self.savePerformance(metrics.SMILESErrorRate.getDjangoModel(), self.current_error)
        self.savePerformance(metrics.MeanDrExActivity.getDjangoModel(), self.current_mean_score)
        if self.best_yet:
            self.savePerformance(metrics.SMILESUniqueRate.getDjangoModel(), self.current_criterion, note=f"Minimum {metrics.SMILESUniqueRate.name}, yet.")
        else:
            self.savePerformance(metrics.SMILESUniqueRate.getDjangoModel(), self.current_criterion)

        if self.original_callback:
            self.original_callback(self)

        self.builder.recordProgress()

    def smiles(self, smiles, score):
        pass

    def performance(self, scores, valids, criterion, best_criterion):
        self.current_criterion = criterion
        self.current_error = 1 - sum(valids) / len(valids)
        self.current_mean_score = scores.mean()

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