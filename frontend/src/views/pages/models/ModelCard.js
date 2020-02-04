import React from "react";
import ModelInfo from './tabs/ModelInfo';
import ModelPerformance from './tabs/ModelPerf';
import ModelPredictions from './tabs/ModelPredictions';
import { ModelCard } from '../../../genui';

class QSARModelCard extends React.Component {

  render() {
    const model =  this.props.model;
    const trainingStrategy =model.trainingStrategy;
    const validationStrategy = model.validationStrategy;

    const trainingParams = [
      {
        name : "Training Set",
        value : model.molset.name
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

    const validationParams = [
      {
        name : "CV-folds",
        value : validationStrategy.cvFolds
      },
      {
        name : "Validation Set Size",
        value : validationStrategy.validSetSize
      }
    ];

    const tabs = [
      {
        title : "Info",
        renderedComponent : () =>
          <ModelInfo
            {...this.props}
            extraTrainingParams={trainingParams}
            extraValidationParams={validationParams}
          />
      },
      {
        title: "Performance"
        , renderedComponent : () =>
          <ModelPerformance
            {...this.props}
          />
      },
      {
        title: "Predictions"
        , renderedComponent : () =>
          <ModelPredictions
            {...this.props}
          />
      }
    ];

    return <ModelCard {...this.props} tabs={tabs}/>
  }
}

export default QSARModelCard;