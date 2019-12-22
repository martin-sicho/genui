"""
serializers

Created by: Martin Sicho
On: 12/22/19, 6:26 PM
"""
from rest_framework import serializers
from rest_framework.schemas.openapi import AutoSchema


class TaskSerializer(serializers.Serializer):
    task_id = serializers.CharField(allow_blank=False)
    status = serializers.CharField(allow_blank=False)

class TasksSerializerFactory:
    class Schema(AutoSchema):
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
        return type('MolSetTasksSerializer', (serializers.Serializer,), {x : serializers.ListField(child=TaskSerializer(required=False), allow_empty=True) for x in field_names})

