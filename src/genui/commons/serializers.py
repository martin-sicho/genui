"""
serializers

Created by: Martin Sicho
On: 12/22/19, 6:26 PM
"""
from rest_framework import serializers


class TaskSerializer(serializers.Serializer):
    task_id = serializers.CharField(allow_blank=False)
    status = serializers.CharField(allow_blank=False)
    result = serializers.JSONField(allow_null=True)
    traceback = serializers.CharField(allow_null=True, max_length=100000)

class TasksSerializerFactory:
    class AutoSchemaMixIn:
        def get_operation(self, path, method):
            ret = super().get_operation(path, method)
            ret['responses']['200']['content']['application/json']['schema'] = {
                'type' : 'object',
                'properties' : {
                        'taskName' : {
                            'type' : 'array',
                            'items' : {
                                'properties' : {
                                    'task_id' : {
                                        'type' : 'string'
                                    },
                                    'status' : {
                                        'type' : 'string'
                                    },
                                    'required': ['task_id', 'status']
                                }
                            }
                        },
                        'required': ['taskName']
                    }
                }
            return ret

    @staticmethod
    def get(field_names):
        return type('TaskSetSerializer', (serializers.Serializer,), {x : serializers.ListField(child=TaskSerializer(required=False), allow_empty=True) for x in field_names})

class TaskProgressInfoSerializer(serializers.Serializer):
    current = serializers.IntegerField()
    total = serializers.IntegerField()
    percent = serializers.FloatField()

class TaskProgressSerializer(serializers.Serializer):
    success = serializers.BooleanField(allow_null=True)
    complete = serializers.BooleanField()
    progress = TaskProgressInfoSerializer()

class GenericModelSerializerMixIn:
    className = serializers.CharField(required=False, allow_blank=True)
    extraArgs = serializers.DictField(required=False, allow_empty=True)

    def to_representation(self, instance):
        ret = super().to_representation(instance)
        ret.update(self.getClassNameRepresentation(instance))
        ret.update(self.getExtraFieldsRepresentation(instance))
        return ret

    def getClassNameRepresentation(self, instance):
        ret = dict()
        model_class = instance.__class__
        if "className" in self.fields.keys():
            ret['className'] = model_class.__name__
        return ret

    def getExtraFieldsRepresentation(self, instance):
        ret = dict()
        if "extraArgs" in self.fields.keys():
            base_fields = set(self.Meta.model._meta.get_fields())
            derived_fields = set(instance.__class__._meta.get_fields())
            extra_fields = []
            for x in derived_fields - base_fields:
                if not x.name.endswith("_ptr"):
                    extra_fields.append(x.name)

            if extra_fields:
                serializer_class = self.getBaseSerializerClass(instance, extra_fields)
                extra_data = serializer_class(instance).data
                ret['extraArgs'] = extra_data
            else:
                ret['extraArgs'] = {}
        return ret

    def getBaseSerializerClass(self, instance, extraFields):
        class GenericModelSerializer(serializers.ModelSerializer):

            class Meta:
                model = instance.__class__
                fields = extraFields
        return GenericModelSerializer