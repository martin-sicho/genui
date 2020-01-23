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
    metrics.forEach(
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
          if (perfMatrix[key].length <= i) {
            retItem[key] = "unavailable";
            retItem.fold = i;
            return
          }
          retItem[key] = perfMatrix[key][i].value;
          retItem["fold"] = perfMatrix[key][i].extraArgs.fold;
        }
      );
      ret.push(retItem);
    }
    ret.sort((a, b) => a.fold <= b.fold ? -1 : 1);
    Object.keys(perfMatrix).forEach(key => {
      const tmp = {fold: "MIN"};
      const arr = perfMatrix[key].map(x => x.value);
      tmp[key] = Math.min(...arr);
      ret.push(tmp);
    });
   Object.keys(perfMatrix).forEach(key => {
      const tmp = {fold: "MAX"};
      const arr = perfMatrix[key].map(x => x.value);
      tmp[key] = Math.max(...arr);
      ret.push(tmp);
    });
    Object.keys(perfMatrix).forEach(key => {
      const tmp = {fold: "AVG"};
      const arr = perfMatrix[key].map(x => x.value);
      tmp[key] = arr.reduce((a,b) => a + b, 0) / arr.length;
      if (Number.isNaN(tmp[key])) {
        tmp[key] = NaN.toString();
      }
      ret.push(tmp);
    });
    Object.keys(perfMatrix).forEach(key => {
      const tmp = {fold: "SD"};
      const arr = perfMatrix[key].map(x => x.value);
      const m = arr.reduce((a,b) => a + b, 0) / arr.length;
      tmp[key] =Math.sqrt(arr.reduce((sq, n) => {
        return sq + Math.pow(n - m, 2);
      }, 0) / (arr.length - 1));
      ret.push(tmp);
    });
    return ret;
  };

  render() {
    const model = this.props.model;
    const validationInfo =  model.validationStrategy;
    const performanceInfo = model.performance;
    const metrics = validationInfo.metrics;
    const metricsNames = metrics.map((metric) => metric.name);

    let validationPerf = this.getPerfMatrix(performanceInfo, "ModelPerformance", metrics);
    validationPerf = Object.keys(validationPerf).map((x) => validationPerf[x].length > 0 ? validationPerf[x][0] : null);
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
        {
          validationPerf[0] !== null ? (
            <Table size="sm">
              <TableDataFromItems
                items={validationPerf}
                dataProps={["value"]}
                rowHeaderProp="metric.name"
              />
            </Table>
          ) : <div>No data.</div>
        }


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