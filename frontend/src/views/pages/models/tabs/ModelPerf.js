import React from "react";
import { Col, Row, Table } from 'reactstrap';
import { TableDataFromItems, TableHeaderFromItems } from '../../../../genui';

class TrainingSummary extends React.Component {

  getPerfValuesForMetric = (performanceInfo, className, metric) => {
    const ret = [];
    performanceInfo.forEach(
      perf => {
        if ((metric ? perf.metric.id === metric.id : true)
          && perf.className === className) {
          ret.push(perf);
        }
      }
    );
    return ret;
  };

  getPerfMatrix = (performanceInfo, className, metrics) => {
    const ret = {};
    metrics.map(
      metric => {
        ret[metric.name] = this.getPerfValuesForMetric(performanceInfo, className, metric);
      });
    return ret;
  };

  parseCVData = (perfMatrix, nFolds) => {
    const ret = [];
    for (let i = 0; i < nFolds; i++) {
      const retItem = {};
      Object.keys(perfMatrix).forEach(
        (key) => {
          retItem[key] = perfMatrix[key][i].value;
          retItem["fold"] = perfMatrix[key][i].extraArgs.fold;
        }
      );
      ret.push(retItem);
    }
    ret.sort((a, b) => a.fold <= b.fold ? -1 : 1);
    return ret;
  };

  render() {
    const model = this.props.model;
    const validationInfo =  model.validationStrategy;
    const performanceInfo = model.performance;
    const metrics = validationInfo.metrics;
    const metricsNames = metrics.map((metric) => metric.name);

    let validationPerf = this.getPerfMatrix(performanceInfo, "ModelPerformance", metrics);
    let cvPerf = this.getPerfMatrix(performanceInfo, "ModelPerformanceCV", metrics);

    return (
      <React.Fragment>
        <h4>
          Training Summary
        </h4>

        <h5>Cross-Validation</h5>
        <Table size="sm" hover striped>
          <TableHeaderFromItems
            items={["Fold"].concat(metricsNames)}
          />
          <TableDataFromItems
            items={this.parseCVData(cvPerf, validationInfo.cvFolds)}
            dataProps={metricsNames}
            rowHeaderProp="fold"
          />
        </Table>

        <h5>Independent Validation Set</h5>
        <Table size="sm">
          <TableDataFromItems
            items={metrics.map((metric) => {
              metric.value = validationPerf[metric.name][0].value; // TODO: check for index out of bounds
              return metric;
            })}
            dataProps={["value"]}
            rowHeaderProp="name"
          />
        </Table>


      </React.Fragment>
    )
  }
}

class ModelPerformance extends React.Component {
  render() {
    const model = this.props.model;

    return (<Row>
      <Col sm="12">
        <TrainingSummary model={model}/>
      </Col>
    </Row>)
  }
}

export default ModelPerformance;