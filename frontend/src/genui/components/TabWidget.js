import React, {Component} from 'react';
import classnames from 'classnames';
import {
  TabContent,
  TabPane,
  Nav,
  NavItem,
  NavLink,
  Card,
} from 'reactstrap';

class TabWidget extends Component {
  constructor(props) {
    super(props);

    this.toggle = this.toggle.bind(this);
    this.state = {
      activeTab: this.props.activeTab
    };
  }

  toggle(tab) {
    if (this.state.activeTab !== tab) {
      this.setState({
        activeTab: tab
      });
    }
  }

  componentDidUpdate(prevProps, prevState, snapshot) {
    if (this.props.activeTab && (prevProps.activeTab !== this.props.activeTab)) {
      this.toggle(this.props.activeTab);
    }
  }

  render() {
    const tabs = this.props.tabs;
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
          <Nav tabs>
            {
              tabs.map(tab => (
                <NavItem key={tab.title}>
                  <NavLink
                    className={classnames({ active: activeTab === tab.title })}
                    onClick={() => { this.toggle(tab.title); }}
                  >
                    {tab.title}
                  </NavLink>
                </NavItem>
              ))
            }
          </Nav>
          <TabContent activeTab={activeTab}>
            {
              tabs.map(tab => {
                const Component = tab.renderedComponent;
                return (
                  <TabPane key={tab.title} tabId={tab.title}>
                    <Component {...this.props}/>
                  </TabPane>
                )
              })
            }
          </TabContent>
        </div>
      </Card>
    )
  }
}

export default TabWidget;