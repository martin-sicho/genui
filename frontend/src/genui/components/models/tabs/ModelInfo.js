import React from "react";
import { Col, Row, Table } from 'reactstrap';
import { TableDataFromItems, TableHeaderFromItems, TaskBadgeGroup, TaskProgressBar, DownloadFile } from '../../../index';

class ModelInfo extends React.Component {

  constructor(props) {
    super(props);

    this.model = this.props.model;
    const trainingStrategy =  this.model.trainingStrategy;
    const validationStrategy = this.model.validationStrategy;

    this.trainingParams = [
      {
        name : "Algorithm",
        value : trainingStrategy.algorithm.name
      },
      {
        name : "Parameters",
        value : trainingStrategy.parameters.map((param) => `${param.parameter.name}=${param.value}`).join(";")
      },
      {
        name : "Mode",
        value : trainingStrategy.mode.name
      },
    ].concat(
      this.props.extraTrainingParams
    );

    this.validationParams = [
      {
        name: "Metrics",
        value: validationStrategy.metrics.map((metric) => `${metric.name}`).join(";")
      }
    ].concat(this.props.extraValidationParams);
  }

  render() {
    const model = this.model;
    const tasks = this.props.tasks;

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
              items={this.trainingParams}
              dataProps={["value"]}
              rowHeaderProp="name"
            />
          </Table>

          <h4>Validation Settings</h4>
          <Table size="sm">
            <TableHeaderFromItems
              items={["Parameter", "Value"]}
            />
            <TableDataFromItems
              items={this.validationParams}
              dataProps={["value"]}
              rowHeaderProp="name"
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