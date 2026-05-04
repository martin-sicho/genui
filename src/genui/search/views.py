from django.shortcuts import render

# Create your views here.

"""
views.py in src/genui/search/

Views of the search package.
"""

from rest_framework.decorators import action
from rest_framework.response import Response
from rest_framework import status
from rest_framework.generics import GenericAPIView
from rest_framework.viewsets import GenericViewSet

from django_rdkit.models import *
from django.contrib.postgres.aggregates import ArrayAgg
from django.db.models.functions import JSONObject

from genui.compounds.models import Molecule, MolSet
from genui.projects.models import Project

from genui.search.serializers import (
    SimilaritySearchRequestSerializer, 
    SimilaritySearchResponseSerializer, 
    SubstructureSearchRequestSerializer,
    SubstructureSearchResponseSerializer,
    SmartsSearchRequestSerializer,
    SmartsSearchResponseSerializer,
    # PropertyFilterSerializer,
    # PropertyFiltersResponseSerializer,
    BaseSearchRequestSerializer,
    BaseSearchResponseSerializer,
    InchiKeySearchRequestSerializer,
    OccurrenceSearchResponseSerializer,
    InchiKeySearchResponseSerializer
    )


class BaseSearch(GenericAPIView):
    queryset = Molecule.objects.all()
    serializer_class = BaseSearchRequestSerializer
    response_serializer_class = BaseSearchResponseSerializer

    def get_molset_ids(self, request, ids):
        return None
    
    def build_queryset(self, molset_ids, params):
        raise NotImplementedError
    
    def serialize_params(self, request):
        params_ser = self.serializer_class(data=request.data)
        params_ser.is_valid(raise_exception=True)
        params = params_ser.validated_data
        return params
    
    def get_searched_count(self, request, molset_ids):
        if molset_ids is None:
            condition = Q(providers__project__owner=request.user)
        else:
            condition = Q(providers__id__in=molset_ids) & Q(providers__project__owner=request.user)

        total_searched = self.queryset.filter(condition).distinct().count()
        return total_searched
    
    def build_response(self,params,hits,total_searched):
        resp = {
            "query": params,
            "hits": hits,
            "total_searched": total_searched,
            "total_returned": len(hits),
        }
        return resp
    
    def do_search(self, request):
        
        params = self.serialize_params(request)
        ids = params["ids"] if "ids" in params else None

        return_value = self.get_molset_ids(request,ids)
        if isinstance(return_value, Response):
            return return_value
        
        molset_ids = return_value
        hits = self.build_queryset(molset_ids,params)
        total_searched = self.get_searched_count(request,molset_ids)
        resp = self.build_response(params,hits,total_searched)

        resp_ser = self.response_serializer_class(resp, context={"request": request})
        return Response(resp_ser.data, status=status.HTTP_200_OK)
    

class SearchMolset(BaseSearch):

    def get_molset_ids(self, request, ids):
        qs = MolSet.objects.filter(project__owner=request.user, pk__in=ids)
        molset_ids = list(qs.values_list("id",flat=True))
        missing = [id for id in ids if id not in molset_ids]
        if missing:
            return Response({"error":f"Some MolSet IDs {missing} not found"}, status=status.HTTP_404_NOT_FOUND)
        return molset_ids
    

class SearchProject(BaseSearch):

    def get_molset_ids(self, request, ids):
        qs = MolSet.objects.filter(project__owner=request.user, project__id__in=ids)
        molset_ids = list(qs.values_list("id",flat=True))
        project_ids = list(Project.objects.filter(id__in=ids).values_list("id",flat=True))
        missing = [id for id in ids if id not in project_ids]
        if missing:
            return Response({"error":f"Some Project IDs {missing} not found"}, status=status.HTTP_404_NOT_FOUND)
        return molset_ids
    

class SimilaritySearch(BaseSearch):

    serializer_class = SimilaritySearchRequestSerializer
    response_serializer_class = SimilaritySearchResponseSerializer

    fingerprints = {
        "maccsFP": MACCS_FP,
        "morganFP": MORGANBV_FP
    }

    sims = {
        "tanimoto":TANIMOTO_SML,
        "dice":DICE_SML
    }

    def build_queryset(self, molset_ids, params):
        smiles = params["canonical"]
        fp_type = params["fp_type"]
        metric = params["metric"]
        threshold = params["threshold"]
        top_n = params["top_n"]

        fp_fn = self.fingerprints[fp_type]
        sim_fn = self.sims[metric]
        value = fp_fn(Value(smiles))

        qs = (
            self.queryset
            .filter(providers__id__in=molset_ids)
            .annotate(similarity=sim_fn(f"entity__{fp_type}", value))
            .order_by("-similarity")
            .filter(similarity__gte=threshold)
            .distinct()
            .annotate(project_ids=ArrayAgg("providers__project__id",distinct=True))
        )

        return list(qs[:top_n])
    

class SubstructureSearch(BaseSearch):
    serializer_class = SubstructureSearchRequestSerializer
    response_serializer_class = SubstructureSearchResponseSerializer

    def build_queryset(self, molset_ids, params):
        smiles = params["canonical"]

        qs = (
            self.queryset
            .filter(Q(providers__id__in=molset_ids) & Q(entity__rdMol__hassubstruct=smiles))
            .order_by(NUMHEAVYATOMS("entity__rdMol"))
            .distinct()
            .annotate(project_ids=ArrayAgg("providers__project__id",distinct=True))
        )

        return list(qs)
    

class SmartsSearch(BaseSearch):
    serializer_class = SmartsSearchRequestSerializer
    response_serializer_class = SmartsSearchResponseSerializer

    def build_queryset(self, molset_ids, params):
        smarts = params["canonical"]

        qs = (
            self.queryset
            .filter(Q(providers__id__in=molset_ids) & Q(entity__rdMol__hassubstruct=QMOL(Value(smarts))))
            .order_by(NUMHEAVYATOMS("entity__rdMol"))
            .distinct()
            .annotate(project_ids=ArrayAgg("providers__project__id",distinct=True))
        )

        return list(qs)
    

class InchiKeySearch(BaseSearch):
    serializer_class = InchiKeySearchRequestSerializer
    response_serializer_class = InchiKeySearchResponseSerializer

    def build_queryset(self, molset_ids, params):
        inchi_key = params["input"]
        qs = (
            self.queryset
            .filter(entity__inchiKey=inchi_key)
        )

        return list(qs)
    

class OccurrenceSearch(BaseSearch):
    queryset = Project.objects.all()
    serializer_class = InchiKeySearchRequestSerializer
    response_serializer_class = OccurrenceSearchResponseSerializer
    
    def build_queryset(self, molset_ids, params):
        inchi_key = params["input"]

        qs = (
            self.queryset
            .filter(molset__molecules__entity__inchiKey=inchi_key)
            .annotate(
                providers=ArrayAgg(
                    JSONObject(
                        id=F("molset__id"),
                        name=F("molset__name"),
                    ),
                    filter=Q(molset__molecules__entity__inchiKey=inchi_key),
                    distinct=True,
                )
            )
        )

        return list(qs)
    
    def get_searched_count(self, request, molset_ids):
        return None
    
    def build_response(self, params, qs, total_searched):
        resp = {
            "query": params,
            "occurrence":qs
        }
        return resp


class SimSearchMolsetView(SearchMolset, SimilaritySearch):

    def post(self, request, *args, **kwargs):
        response = self.do_search(request)
        return response
        
    
class SubsSearchMolsetView(SearchMolset, SubstructureSearch):

    def post(self, request, *args, **kwargs):
        response = self.do_search(request)
        return response

    
class SmartsSearchMolsetView(SearchMolset, SmartsSearch):

    def post(self, request, *args, **kwargs):
        response = self.do_search(request)
        return response


class SimSearchProjectView(SearchProject, SimilaritySearch):

    def post(self, request, *args, **kwargs):
        response = self.do_search(request)
        return response
    

class SubsSearchProjectView(SearchProject, SubstructureSearch):

    def post(self, request, *args, **kwargs):
        response = self.do_search(request)
        return response


class SmartsSearchProjectView(SearchProject, SmartsSearch):

    def post(self, request, *args, **kwargs):
        response = self.do_search(request)
        return response
    

class InchiKeyOccurrenceSearchView(OccurrenceSearch):

    def post(self, request, *args, **kwargs):
        response = self.do_search(request)
        return response
    

class InchiKeyProjectsSearchView(InchiKeySearch):

    def post(self, request, *args, **kwargs):
        response = self.do_search(request)
        return response
    
    
# class PropertyFilterView(SearchMolset):

#     serializer_class = PropertyFilterSerializer
#     response_serializer_class = PropertyFiltersResponseSerializer

#     def build_queryset(self, molset_ids, params):
#         property = params["property"]
#         relation = params["relation"]
#         value = params["value"]

#         qs = (
#             self.get_queryset()
#             .select_related("entity")
#             .prefetch_related("activities")
#             .filter(Q(providers__id__in=molset_ids) & Q(**{f"entity__rdMol__{property}__{relation}":value}))
#             .distinct()
#         )

#         return qs
    
#     def post(self, request, *args, **kwargs):
#         response = self.do_search(request)
#         return response
    