import React from "react"
import { ModelCard, ModelInfoTab, ModelPerformanceTab } from '../../../../genui';

function DrugExPerformanceTab (props) {
  return <div>Performance here...</div>
}

export class DrugExNetCard extends React.Component {

  render() {
    const model =  this.props.model;
    // const trainingStrategy = model.trainingStrategy;
    const validationStrategy = model.validationStrategy;

    const trainingParams = [
      {
        name : "Training Set",
        value : model.molset ? model.molset.name : "---"
      },
      {
        name : "Parent",
        value : model.parent ? model.parent.name : "---"
      },
    ];

    const validationParams = validationStrategy ? [
      {
        name : "Validation Set Size",
        value : validationStrategy.validSetSize
      }
    ] : [];

    const tabs = [
      {
        title : "Info",
        renderedComponent : (props) =>
          <ModelInfoTab
            {...props}
            extraTrainingParams={trainingParams}
            extraValidationParams={validationParams}
          />
      },
      {
        title: "Performance"
        , renderedComponent : (props) =>
          <ModelPerformanceTab
            {...props}
            performanceData={[]} // TODO: fetch performance from server (do this in the tab itself using a url)
            component={DrugExPerformanceTab}
          />
      }
    ];

    return <ModelCard {...this.props} tabs={tabs}/>
  }
}

export class DrugExAgentCard extends React.Component {

  render() {
    const model =  this.props.model;

    const trainingParams = [
      {
        name : "Environment",
        value : model.environment.name
      },
      {
        name : "Exploitation Network",
        value : model.exploitationNet.name
      },
      {
        name : "Exploration Network",
        value : model.explorationNet.name
      },
    ];

    const tabs = [
      {
        title : "Info",
        renderedComponent : () =>
          <ModelInfoTab
            {...this.props}
            extraTrainingParams={trainingParams}
          />
      },
      {
        title: "Performance"
        , renderedComponent : () =>
          <ModelPerformanceTab
            {...this.props}
            performanceData={[]} // TODO: fetch performance from server (do this in the tab itself using a url)
            component={DrugExPerformanceTab}
          />
      }
    ];

    return <ModelCard {...this.props} tabs={tabs}/>
  }
}