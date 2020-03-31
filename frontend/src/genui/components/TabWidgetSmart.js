import React from 'react';
import { Tab, Tabs } from 'react-bootstrap';
import classnames from 'classnames';
import { Card } from 'reactstrap';

class TabWidgetSmart extends React.Component {
  constructor(props) {
    super(props);

    const tabs = this.props.tabs;
    const activeTab = this.props.activeTab ? this.props.activeTab : 'Info';
    this.toggle = this.toggle.bind(this);
    this.state = {
      activeTab: activeTab,
      tabs: tabs
    }
  }

  toggle(tab) {
    if (this.state.activeTab !== tab) {
      this.setState({
        activeTab: tab
      });
    }
  }
  render() {
    const tabs = this.state.tabs;
    let activeTab = tabs.find(tab => tab.title === this.state.activeTab);
    if (!activeTab) {
      activeTab = tabs[0]
    }
    if (!activeTab) {
      throw new Error("No valid active tab found. Make sure to specify it in a prop or in the state.");
    }
    activeTab = activeTab.title;

    return (
      <Card body className="stretch-to-container unDraggable">
        <div className="full-bleed">
          <Tabs
            defaultActiveKey={tabs[0].title}
            activeKey={activeTab}
            onSelect={k => this.toggle(k)}
            id="noanim-tab-example"
            unmountOnExit={false}
            mountOnEnter={false}
          >
            {
              tabs.map(tab => {
                const Component = tab.renderedComponent;
                return (
                  <Tab
                    eventKey={tab.title}
                    key={tab.title}
                    title={tab.title}
                    className={classnames({ active: activeTab === tab.title })}
                  >
                    <Component {...this.props}/>
                  </Tab>
                )
              })
            };
            }
          </Tabs>
        </div>
      </Card>
    )
  }
}

export default TabWidgetSmart;