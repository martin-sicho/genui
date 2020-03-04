import React from "react";
import { Col, Row} from 'reactstrap';

class ModelPerformance extends React.Component {

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

  render() {
    const SummaryComponent = this.props.component;

    return (<Row>
      <Col sm="12">
        <SummaryComponent {...this.props} getPerfMatrix={this.getPerfMatrix} getPerfValuesForMetric={this.getPerfValuesForMetric}/>
      </Col>
    </Row>)
  }
}

export default ModelPerformance;