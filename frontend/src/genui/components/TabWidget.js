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
      activeTab: this.props.activeTab ? this.props.activeTab : 'Info'
    };
  }

  toggle(tab) {
    if (this.state.activeTab !== tab) {
      this.setState({
        activeTab: tab
      });
    }
  }
  render() {
    const tabs = this.props.tabs;

    return (
      <Card body className="stretch-to-container unDraggable">
        <div className="full-bleed">
          <Nav tabs>
            {
              tabs.map(tab => (
                <NavItem key={tab.title}>
                  <NavLink
                    className={classnames({ active: this.state.activeTab === tab.title })}
                    onClick={() => { this.toggle(tab.title); }}
                  >
                    {tab.title}
                  </NavLink>
                </NavItem>
              ))
            }
          </Nav>
          <TabContent activeTab={this.state.activeTab}>
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