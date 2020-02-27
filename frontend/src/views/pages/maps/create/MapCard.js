import React from 'react';
import { ModelCard, ModelInfoTab } from '../../../../genui';

export default function MapCard(props) {
  const model =  props.model;
  const trainingStrategy = model.trainingStrategy;

  const trainingParams = [
    {
      name : "Training Set",
      value : model.molsets.map(molset => `${molset.name}`).join(";")
    },
    {
      name : "Descriptor Sets",
      value : trainingStrategy.descriptors.map((desc) => `${desc.name}`).join(";")
    }
  ];

  const tabs = [
    {
      title : "Info",
      renderedComponent : () =>
        <ModelInfoTab
          {...props}
          extraTrainingParams={trainingParams}
        />
    }
  ];

  return <ModelCard {...props} tabs={tabs}/>
}