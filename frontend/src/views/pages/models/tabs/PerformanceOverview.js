import React from 'react';
import { Table } from 'reactstrap';
import { LiveObject, TableDataFromItems, TableHeaderFromItems } from '../../../../genui';

class QSARPerformance extends React.Component {

  parseCVData = (perfMatrix, nFolds) => {
    const ret = [];
    for (let i = 0; i < nFolds; i++) {
      const retItem = {};
      Object.keys(perfMatrix).forEach(
        (key) => {
          if (perfMatrix[key].length <= i) {
            retItem[key] = 'unavailable';
            retItem.fold = i;
            return;
          }
          retItem[key] = perfMatrix[key][i].value;
          retItem['fold'] = perfMatrix[key][i].extraArgs.fold;
        },
      );
      ret.push(retItem);
    }
    ret.sort((a, b) => a.fold <= b.fold ? -1 : 1);
    Object.keys(perfMatrix).forEach(key => {
      const tmp = { fold: 'MIN' };
      const arr = perfMatrix[key].map(x => x.value);
      tmp[key] = Math.min(...arr);
      ret.push(tmp);
    });
    Object.keys(perfMatrix).forEach(key => {
      const tmp = { fold: 'MAX' };
      const arr = perfMatrix[key].map(x => x.value);
      tmp[key] = Math.max(...arr);
      ret.push(tmp);
    });
    Object.keys(perfMatrix).forEach(key => {
      const tmp = { fold: 'AVG' };
      const arr = perfMatrix[key].map(x => x.value);
      tmp[key] = arr.reduce((a, b) => a + b, 0) / arr.length;
      if (Number.isNaN(tmp[key])) {
        tmp[key] = NaN.toString();
      }
      ret.push(tmp);
    });
    Object.keys(perfMatrix).forEach(key => {
      const tmp = { fold: 'SD' };
      const arr = perfMatrix[key].map(x => x.value);
      const m = arr.reduce((a, b) => a + b, 0) / arr.length;
      tmp[key] = Math.sqrt(arr.reduce((sq, n) => {
        return sq + Math.pow(n - m, 2);
      }, 0) / (arr.length - 1));
      ret.push(tmp);
    });
    return ret;
  };

  render() {
    const model = this.props.model;
    const validationStratInfo = model.validationStrategy;
    if (!validationStratInfo) {
      return <p>No performance data for this model is available.</p>
    }

    const performanceInfo = model.performance;
    const metrics = validationStratInfo.metrics;
    const metricsNames = metrics.map((metric) => metric.name);

    let validSetPerf = this.props.getPerfMatrix(performanceInfo, 'ModelPerformance', metrics);
    validSetPerf = Object.keys(validSetPerf).map((x) => validSetPerf[x].length > 0 ? validSetPerf[x][0] : null);
    let cvPerf = this.props.getPerfMatrix(performanceInfo, 'ModelPerformanceCV', metrics);
    return (
      <React.Fragment>
        <h4>
          Training Summary
        </h4>

        <h5>Cross-Validation</h5>
        <Table size="sm" hover striped>
          <TableHeaderFromItems
            items={['Fold'].concat(metricsNames)}
          />
          <TableDataFromItems
            items={this.parseCVData(cvPerf, validationStratInfo.cvFolds)}
            dataProps={metricsNames}
            rowHeaderProp="fold"
          />
        </Table>

        <h5>Independent Validation Set</h5>
        {
          validSetPerf[0] !== null ? (
            <Table size="sm">
              <TableDataFromItems
                items={validSetPerf}
                dataProps={['value']}
                rowHeaderProp="metric.name"
              />
            </Table>
          ) : <div>No data.</div>
        }
      </React.Fragment>
    );
  }
}

function QSARPerformanceLive(props) {
  return (
    <LiveObject {...props} url={props.modelUrl}>
      {
        (model) => (
          <QSARPerformance {...props} model={model}/>
        )
      }
    </LiveObject>
  )
}

export function QSARPerformanceOverview(props) {
  if (props.tasksRunning) {
    return <QSARPerformanceLive {...props}/>
  } else {
    return <QSARPerformance {...props}/>
  }
}

export default QSARPerformanceOverview;