import React from "react";
import ModelPerformance from './tabs/ModelPerf';
import ModelPredictions from './tabs/ModelPredictions';
import { ModelCard, ModelInfoTab, TableDataFromItems, TableHeaderFromItems } from '../../../genui';
import { Table } from 'reactstrap';

class QSARPerformanceSummary extends React.Component {

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
    const performanceInfo = this.props.performanceData;
    const metrics = validationInfo.metrics;
    const metricsNames = metrics.map((metric) => metric.name);

    let validationPerf = this.props.getPerfMatrix(performanceInfo, "ModelPerformance", metrics);
    validationPerf = Object.keys(validationPerf).map((x) => validationPerf[x].length > 0 ? validationPerf[x][0] : null);
    let cvPerf = this.props.getPerfMatrix(performanceInfo, "ModelPerformanceCV", metrics);

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
          <ModelInfoTab
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
            performanceData={model.performance}
            component={QSARPerformanceSummary}
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