import React from 'react';
import { ModelCard, ModelInfoTab, ModelPerformanceTab, ModelPreds } from '../../../genui';
import QSARPerformanceOverview from './tabs/PerformanceOverview';

class QSARModelCard extends React.Component {

  render() {
    const model =  this.props.model;
    const trainingStrategy = model.trainingStrategy;

    const modelData = [
      {
        name : "Compound Set",
        value : model.molset ? model.molset.name : ""
      },
      {
        name : "Predictions Activity Type",
        value : model.predictionsType ? model.predictionsType.value : ""
      },
      {
        name : "Predictions Activity Units",
        value : model.predictionsUnits ? model.predictionsUnits.value : ""
      },
    ];

    const trainingParams = [
      {
        name : "Activity Set",
        value : trainingStrategy.activitySet ? trainingStrategy.activitySet.name : ""
      },
      {
        name : "Activity Type",
        value : trainingStrategy.activityType ? trainingStrategy.activityType.value : ""
      },
      {
        name : "Activity Threshold",
        value : trainingStrategy.activityThreshold
      },
      {
        name : "Descriptor Sets",
        value : trainingStrategy.descriptors.map((desc) => `${desc.name}`).join(";")
      }
    ];
    if (trainingStrategy.modelledActivityType) {
      trainingParams.push({
        name : "Modelled Activity Type",
        value : trainingStrategy.modelledActivityType.value
      });
      trainingParams.push({
        name : "Modelled Activity Units",
        value : trainingStrategy.modelledActivityUnits ? trainingStrategy.modelledActivityUnits.value : 'No dimension.'
      })
    }

    const validationStrategy = model.validationStrategy;
    const validationParams = [];
    if (validationStrategy) {
      validationParams.push({
        name : "CV-folds",
        value : validationStrategy.cvFolds
      });
      validationParams.push({
          name : "Validation Set Size",
          value : validationStrategy.validSetSize
      });
    }

    const tabs = [
      {
        title : "Info",
        renderedComponent : (props) =>
          <ModelInfoTab
            {...props}
            extraTrainingParams={trainingParams}
            extraValidationParams={validationParams}
            modelData={modelData}
          />
      },
      {
        title: "Performance"
        , renderedComponent : (props) =>
          <ModelPerformanceTab
            {...props}
            component={QSARPerformanceOverview}
          />
      },
      {
        title: "Predictions"
        , renderedComponent : ModelPreds
      }
    ];

    return <ModelCard {...this.props} tabs={tabs}/>
  }
}

export default QSARModelCard;