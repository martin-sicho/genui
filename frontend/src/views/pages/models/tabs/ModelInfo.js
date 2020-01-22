import React from "react";
import { Col, Row, Table } from 'reactstrap';
import { TableDataFromItems, TableHeaderFromItems, TaskBadgeGroup, TaskProgressBar, DownloadFile } from '../../../../genui';

class ModelInfo extends React.Component {

  getTrainingInfo = (trainingStrategy) => {
    const ret = [];
    ret.push({
      paramName : "Training Set",
      value : trainingStrategy.molset.name
    });
    ret.push({
      paramName : "Algorithm",
      value : trainingStrategy.algorithm.name
    });
    ret.push({
      paramName : "Alg. Parameters",
      value : trainingStrategy.parameters.map((param) => `${param.parameter.name}=${param.value}`).join(";")
    });
    ret.push({
      paramName : "Mode",
      value : trainingStrategy.mode.name
    });
    if(trainingStrategy.mode.name === "classification") {
      ret.push({
        paramName : "Activity Threshold",
        value : trainingStrategy.activityThreshold
      })
    }
    ret.push({
      paramName : "Descriptor Sets",
      value : trainingStrategy.descriptors.map((desc) => `${desc.name}`).join(";")
    });
    return ret;
  };

  getValidationInfo = (validationStrategy) => {
    const ret = [];
    ret.push({
      paramName : "CV-folds",
      value : validationStrategy.cvFolds
    });
    ret.push({
      paramName : "Validation Set Size",
      value : validationStrategy.validSetSize
    });
    return ret;
  };

  render() {
    const model = this.props.model;
    const tasks = this.props.tasks;
    const trainingInfo =  model.trainingStrategy;
    trainingInfo.molset = model.molset;
    const validationInfo = model.validationStrategy;

    return (
      (<Row>
        <Col sm="12">
          {
            model.description ? (
              <React.Fragment>
                <h4>Description</h4>
                <p>{model.description}</p>
              </React.Fragment>
            ) : null
          }

          <h4>Training Settings</h4>
          <Table size="sm">
            <TableHeaderFromItems
            items={["Parameter", "Value"]}
            />
            <TableDataFromItems
              items={this.getTrainingInfo(trainingInfo)}
              dataProps={["value"]}
              rowHeaderProp="paramName"
            />
          </Table>

          <h4>Validation Settings</h4>
          <Table size="sm">
            <TableHeaderFromItems
              items={["Parameter", "Value"]}
            />
            <TableDataFromItems
              items={this.getValidationInfo(validationInfo)}
              dataProps={["value"]}
              rowHeaderProp="paramName"
            />
          </Table>

          {
            model.modelFile ? (
              <React.Fragment>
                <h4>Download Model</h4>
                <DownloadFile
                  file={model.modelFile}
                  name={model.modelFile.split("-").slice(-1)[0]}
                />
              </React.Fragment>
            ) : null
          }

          <h4>
            Tasks <TaskBadgeGroup tasks={tasks}/>
          </h4>
          <TaskProgressBar
            progressURL={this.props.apiUrls.celeryProgress}
            tasks={tasks.running}
          />
        </Col>
      </Row>)
    )
  }
}

export default ModelInfo;