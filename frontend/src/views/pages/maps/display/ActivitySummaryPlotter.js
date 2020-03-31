import { groupBy, GroupedViolinPlot } from '../../../../genui';
import React from 'react';

export default function ActivitySummaryPlotter(props) {
  const molsets = props.molsets;
  const actsets = props.activitySets;

  const plotTraces = {};
  const bySource = groupBy(props.activities, 'source');
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
        color: props.molsetsToColor[molset.id],
      }
    }

    const values = group.map(item => item.value);
    const groups = group.map(item => actset.name);
    const customdata = group.map(item => props.selectedMols.find(mol => mol.id === item.molecule));
    plotTraces[molset.id].x = plotTraces[molset.id].x.concat(groups);
    plotTraces[molset.id].y = plotTraces[molset.id].y.concat(values);
    plotTraces[molset.id].customdata = plotTraces[molset.id].customdata.concat(customdata);
  });

  return (
    <GroupedViolinPlot
      {...props}
      title={`Distributions of ${props.type.value}`}
      traces={plotTraces}
      onHover={(eventData) => eventData ? props.onMolHover(eventData.points[0].customdata) : null}
      onSelect={(eventData) => eventData ? props.onMolsSelect(eventData.points.map(point => point.customdata)) : null}
    />
  )
}