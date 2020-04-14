import React from 'react';
import { Col, Row } from 'reactstrap';
import { ActivitiesAggregator, groupBy, TabWidget } from '../../../../genui';
import ActivitySummary from './ActivitySummary';

class SelectedActivitiesPage extends React.Component {

  constructor(props) {
    super(props);

    this.state = {
      currentSelectionRev: props.selectedMolsRevision,
      aggregator: null
    }
  }

  getAggregator = (mols) => {
    return (props) => <ActivitiesAggregator {...props} mols={mols}/>
  };

  componentDidUpdate(prevProps, prevState, snapshot) {
    if (this.props.selectedMols.length > 0 && (this.props.selectedMolsRevision !== this.state.currentSelectionRev)) {
      this.setState({
        currentSelectionRev: this.props.selectedMolsRevision,
        aggregator: this.getAggregator(this.props.selectedMols),
      })
    }
  }

  render() {
    return (
      <React.Fragment>
        <h1>Selected Compounds: Activities Summary</h1>
        <hr/>

        <Row>
          <Col sm={12}>
            {
              this.state.aggregator ? (
                <this.state.aggregator
                  {...this.props}
                >
                  {
                    (activities) => {
                      if (activities) {
                        const groupedActivities = groupBy(activities, 'type.id');
                        // console.log(groupedActivities);

                        const tabs = groupedActivities.map(group => ({
                          title: group[0].type.value,
                          renderedComponent: (props) => (
                            <ActivitySummary
                              {...props}
                              type={group[0].type}
                              activities={group}
                            />
                          )
                        }));

                        return (
                          <TabWidget
                            {...this.props}
                            tabs={tabs}
                          />
                        )
                      } else {
                        return <div>Loading...</div>;
                      }
                    }
                  }
                </this.state.aggregator>
              ) : <p>Select compounds in the map to see details.</p>
            }
          </Col>
        </Row>
      </React.Fragment>
    )
  }
}

export default SelectedActivitiesPage;