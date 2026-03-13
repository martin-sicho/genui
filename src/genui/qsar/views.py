import traceback
from django.conf import settings
from drf_yasg.utils import swagger_auto_schema
from rest_framework import viewsets, mixins, status
from rest_framework.decorators import action, api_view, permission_classes
from rest_framework.permissions import AllowAny
from rest_framework.response import Response
from rest_framework import serializers as rest_serializers

from genui.models import models as genui_models
from genui.models.serializers import HyperparameterOptimizationStrategySerializer
from genui.models.views import ModelViewSet, AlgorithmViewSet, PredictMixIn
from genui.models.genuimodels.bases import Algorithm
from genui.qsar.genuimodels.builders import BasicQSARModelBuilder
from genui.qsar.genuimodels.bases import EmbeddingBuilderMixIn, EmbeddingCalculator
from genui.qsar.genuimodels import embeddings
from . import models
from . import serializers
from .tasks import buildQSARModel, predictWithQSARModel
from genui.utils.extensions.tasks.utils import runTask
from genui import celery_app
from genui.utils.inspection import get_default_params, get_default_params_django, sklearn_regressors, \
    sklearn_classifiers, SKLEARN_MODELS, SKLEARN_MODELS_PARAMS, getSubclassesFromModule, METRICS, DATA_SPLITS, \
    SCAFFOLDS, CLUSTERING


class QSARModelViewSet(PredictMixIn, ModelViewSet):
    queryset = models.QSARModel.objects.order_by('-created')
    serializer_class = serializers.QSARModelSerializer
    init_serializer_class = serializers.QSARModelInitSerializer
    builder_class = BasicQSARModelBuilder
    build_task = buildQSARModel
    predict_task = predictWithQSARModel

    @swagger_auto_schema(
        methods=['GET']
        , responses={
            200: serializers.ModelActivitySetSerializer(many=True),
        }
    )
    @swagger_auto_schema(
        methods=['POST']
        , responses={
            201: serializers.ModelActivitySetSerializer(many=False),
        }
        , request_body=serializers.ModelActivitySetSerializer(many=False)
    )
    @swagger_auto_schema(
        methods=['DELETE']
        , responses={
            204: "",
        }
    )
    @action(detail=True, methods=['get', 'post', 'delete'])
    def predictions(self, request, pk=None):
        # FIXME: some of this should be moved to the PredictMixIn for reuse
        try:
            instance = self.get_queryset().get(pk=pk)
        except models.QSARModel.DoesNotExist:
            return Response({"error": f"Model not found: {pk}"}, status=status.HTTP_404_NOT_FOUND)

        if request.method == 'GET':
            predictions = instance.predictions.all()
            serializer = serializers.ModelActivitySetSerializer(predictions, many=True)
            return Response(serializer.data)

        elif request.method == 'POST':
            request.data['project'] = instance.project.id
            request.data['model'] = instance.id
            serializer = serializers.ModelActivitySetSerializer(data=request.data, many=False)
            if serializer.is_valid():
                created = serializer.create(serializer.validated_data)

                task = None
                try:
                    task, task_id = runTask(
                        self.get_predict_task(),
                        instance=instance,
                        eager=hasattr(settings, 'CELERY_TASK_ALWAYS_EAGER') and settings.CELERY_TASK_ALWAYS_EAGER,
                        args=(
                            created.pk,
                            self.get_builder_class()
                        )
                    )
                    ret = serializers.ModelActivitySetSerializer(created, many=False)
                    ret.data["taskID"] = task_id
                    return Response(ret.data, status=status.HTTP_201_CREATED)
                except Exception as exp:
                    traceback.print_exc()
                    if task and task.id:
                        celery_app.control.revoke(task_id=task.id, terminate=True)
                    created.delete()
                    return Response({"error": repr(exp)}, status=status.HTTP_500_INTERNAL_SERVER_ERROR)

            return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)

        elif request.method == 'DELETE':
            instance.predictions.all().delete()
            return Response(status=status.HTTP_204_NO_CONTENT)
        else:
            return Response({"error": "Method not allowed"}, status=status.HTTP_405_METHOD_NOT_ALLOWED)


class EmbeddingGroupsViewSet(
    mixins.ListModelMixin,
    mixins.RetrieveModelMixin,
    mixins.CreateModelMixin,
    viewsets.GenericViewSet
):
    queryset = models.EmbeddingCalculator.objects.all()
    serializer_class = serializers.EmbeddingCalculatorSerializer

    @swagger_auto_schema(
        methods=['GET']
        , responses={
            200: "List of all available embedding calculator names",
        }
    )
    @action(detail=False, methods=['get'], url_path='list')
    def list_all(self, request):
        calculators = getSubclassesFromModule(EmbeddingCalculator, embeddings)
        calculators = [c.__name__ for c in calculators if not c.abstract]
        return Response(calculators)

    @action(detail=False, methods=['get'], url_path='(?P<name>[^/.]+)/arguments')
    def params_by_name(self, request, name=None):
        class_ = EmbeddingBuilderMixIn.findEmbeddingClass(name)
        if not class_:
            return Response({"error": f"Embedding class not found for: {name}"},
                            status=status.HTTP_404_NOT_FOUND)

        calculator = class_(None)
        data = calculator.get_default_parameters()
        return Response(data)


class QSARAlgorithmViewSet(AlgorithmViewSet):

    def get_queryset(self):
        current = super().get_queryset()
        return current.filter(validModes__name__in=(Algorithm.CLASSIFICATION, Algorithm.REGRESSION)).distinct('id')


@api_view(['GET'])
@permission_classes([AllowAny])
def list_metrics(request):
    return Response(METRICS)


@api_view(['GET'])
@permission_classes([AllowAny])
def list_metrics_by_mode(request, mode):
    if mode == Algorithm.REGRESSION:
        metrics = METRICS[Algorithm.REGRESSION]
    elif mode == Algorithm.CLASSIFICATION:
        metrics = METRICS[Algorithm.CLASSIFICATION]
    else:
        return Response({"error": f"Mode not found: {mode}"}, status=status.HTTP_404_NOT_FOUND)
    serializer = rest_serializers.ListSerializer(data=metrics, child=rest_serializers.CharField())
    serializer.is_valid(raise_exception=True)
    return Response(serializer.data)


@api_view(['GET'])
@permission_classes([AllowAny])
def list_data_splits(request):
    serializer = rest_serializers.ListSerializer(data=DATA_SPLITS, child=rest_serializers.CharField())
    serializer.is_valid(raise_exception=True)
    return Response(serializer.data)


@api_view(['GET'])
@permission_classes([AllowAny])
def list_aggregation_functions(request):
    serializer = rest_serializers.ListSerializer(data=["max", "min", "mean", "sum", "median", "std", "var"],
                                                 child=rest_serializers.CharField())
    serializer.is_valid(raise_exception=True)
    return Response(serializer.data)


@api_view(['GET'])
@permission_classes([AllowAny])
def list_scaffolds(request):
    serializer = rest_serializers.ListSerializer(data=SCAFFOLDS, child=rest_serializers.CharField())
    serializer.is_valid(raise_exception=True)
    return Response(serializer.data)


@ api_view(['GET'])
@permission_classes([AllowAny])
def list_clustering(request):
    serializer = rest_serializers.ListSerializer(data=CLUSTERING, child=rest_serializers.CharField())
    serializer.is_valid(raise_exception=True)
    return Response(serializer.data)


@api_view(['GET'])
@permission_classes([AllowAny])
def get_data_splits_params(request, data_split):
    module_name, func_name = data_split.rsplit(".", 1)
    args = get_default_params(func_name, module_name)
    unwanted_args = ["dataset", "weights", "custom"]
    args = {k: v for k, v in args.items() if not any(sub in k for sub in unwanted_args)}
    return Response(args)


class QSPRPredSklearnModelViewSet(viewsets.ViewSet):
    queryset = models.QSPRPredSklearnModel.objects.all()
    serializer_class = serializers.QSPRPredSklearnModelSerializer

    def list(self, request):
        data = [name for name in SKLEARN_MODELS]
        serializer = rest_serializers.ListSerializer(data=data, child=rest_serializers.CharField())
        serializer.is_valid(raise_exception=True)
        return Response(serializer.data)

    def mode_model(self, request, mode):
        if mode == "regression":
            model_list = [name for name in sklearn_regressors]
        elif mode == "classification":
            model_list = [name for name in sklearn_classifiers]
        else:
            return Response({"error": f"Mode not found: {mode}"}, status=status.HTTP_404_NOT_FOUND)
        serializer = rest_serializers.ListSerializer(data=model_list, child=rest_serializers.CharField())
        serializer.is_valid(raise_exception=True)
        return Response(serializer.data)

    def retrieve(self, request, algorithm):
        data = []
        if algorithm in sklearn_regressors:
            cls = sklearn_regressors[algorithm]
            data.append({'name': algorithm, 'type': 'Regressor', 'params': get_default_params(None, cls)})
        elif algorithm in sklearn_classifiers:
            cls = sklearn_classifiers[algorithm]
            data.append({'name': algorithm, 'type': 'Classifier', 'params': get_default_params(None, cls)})
        else:
            return Response({"error": f"Algorithm not found: {algorithm}"}, status=status.HTTP_404_NOT_FOUND)
        serializer = serializers.QSPRPredSklearnModelSerializer(data, many=True)
        return Response(serializer.data)

    def get_type(self, request, algorithm):
        if algorithm in sklearn_regressors:
            return Response("Regressor")
        elif algorithm in sklearn_classifiers:
            return Response("Classifier")
        else:
            return Response({"error": f"Algorithm not found: {algorithm}"}, status=status.HTTP_404_NOT_FOUND)

    def get_params(self, request, algorithm):
        if algorithm in SKLEARN_MODELS:
            data = SKLEARN_MODELS_PARAMS[algorithm]
        else:
            return Response({"error": f"Algorithm not found: {algorithm}"}, status=status.HTTP_404_NOT_FOUND)
        return Response(data)


class QSARHyperParameterOptimizationViewSet(viewsets.ViewSet):
    queryset = genui_models.HyperparameterOptimizationStrategy.objects.all()
    serializer_class = HyperparameterOptimizationStrategySerializer

    @swagger_auto_schema(
        methods=['GET']
        , responses={
            200: "List of all available hyperparameter optimization strategies",
        }
    )
    @action(detail=False, methods=['get'], url_path='list')
    def list_all(self, request):
        from django.contrib.contenttypes.models import ContentType
        from genui.models.models import HyperparameterOptimizationStrategy

        hyperparam_strats = []
        for ct in ContentType.objects.filter(app_label__in=['models']):
            model_class = ct.model_class()
            if model_class and issubclass(model_class,
                                          HyperparameterOptimizationStrategy) and model_class != HyperparameterOptimizationStrategy:
                hyperparam_strats.append(model_class.__name__)

        return Response(hyperparam_strats)

    @swagger_auto_schema(
        methods=['GET']
        , responses={
            200: "Default parameters for the hyperparameter optimization strategy",
        }
    )
    @action(detail=False, methods=['get'], url_path='(?P<name>[^/.]+)/params')
    def params_by_name(self, request, name=None):
        from django.contrib.contenttypes.models import ContentType
        from genui.models.models import HyperparameterOptimizationStrategy

        model_class = None
        for ct in ContentType.objects.filter(app_label__in=['models']):
            cls = ct.model_class()
            if cls and issubclass(cls, HyperparameterOptimizationStrategy) and cls.__name__ == name:
                model_class = cls
                break

        if not model_class:
            return Response({"error": f"Hyperparameter optimization strategy not found: {name}"},
                            status=status.HTTP_404_NOT_FOUND)

        try:
            module_name = model_class.__module__
            class_name = model_class.__name__
            params = get_default_params_django(class_name, module_name)
            return Response(params)
        except Exception as e:
            return Response({"error": f"Failed to get parameters: {str(e)}"},
                            status=status.HTTP_500_INTERNAL_SERVER_ERROR)
