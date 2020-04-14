import { groupBy, GroupedViolinPlot } from '../../../../genui';
import React from 'react';

class ActivitySummaryPlotter extends React.Component {

  constructor(props) {
    super(props);

    this.state = {
      traces: this.initTraces(),
      revision: 0
    }
  }

  initTraces = () => {
    const molsets = this.props.molsets;
    const actsets = this.props.activitySets;
    const activities = this.props.activities;
    const mols = this.props.mols;

    const plotTraces = {};
    const bySource = groupBy(activities, 'source');
    bySource.forEach(group => {
      const actset = actsets[group[0].source];
      const molset = molsets.find(item => item.id === actset.molecules);

      if (!plotTraces.hasOwnProperty(molset.id)) {
        plotTraces[molset.id] = {};
        plotTraces[molset.id].name = molset.name;
        plotTraces[molset.id].x = [];
        plotTraces[molset.id].y = [];
        plotTraces[molset.id].customdata = [];
        plotTraces[molset.id].marker = {
          color: this.props.molsetsToColor[molset.id],
        }
      }

      const values = group.map(item => item.value);
      const groups = group.map(item => actset.name);
      const customdata = group.map(item => mols.find(mol => mol.id === item.molecule));
      plotTraces[molset.id].x = plotTraces[molset.id].x.concat(groups);
      plotTraces[molset.id].y = plotTraces[molset.id].y.concat(values);
      plotTraces[molset.id].customdata = plotTraces[molset.id].customdata.concat(customdata);
    });

    return plotTraces;
  };

  componentDidUpdate(prevProps, prevState, snapshot) {
    if (prevProps.selectedMolsRevision !== this.props.selectedMolsRevision) {
      this.setState(prevState => ({
        traces: this.initTraces(),
        revision: prevState.revision + 1
      }))
    }
  }

  render() {
    return (
      <GroupedViolinPlot
        title={`Distributions of ${this.props.type.value}`}
        traces={this.state.traces}
        tracesRev={this.state.revision}
        onHover={(eventData) => eventData ? this.props.onMolHover(eventData.points[0].customdata) : null}
        onSelect={(eventData) => {
          if (eventData) {
            this.props.onMolsSelect(eventData.points.map(point => point.customdata));
            if (this.props.setSelectedMolsOverviewRevision) {
              this.props.setSelectedMolsOverviewRevision(this.props.selectedMolsOverviewRevision + 1);
            }
          }
        }
        }
      />
    )
  }
}

export default ActivitySummaryPlotter;